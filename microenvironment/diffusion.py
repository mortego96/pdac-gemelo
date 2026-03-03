"""
Microambiente tumoral PDAC v5.
9 campos de difusión + lactato, pH, colágeno, TNFα.
Basado en datos de PhysiCell, modelos PDAC publicados.
"""

import numpy as np
from scipy.ndimage import uniform_filter


class Microenvironment:
    """
    Gestiona los campos espaciales del microambiente tumoral.
    Cada campo es un array 2D (width × height), valores [0, 1].
    """

    def __init__(self, width: int, height: int):
        self.width = width
        self.height = height

        # ═══ CAMPOS PRIMARIOS ═══
        self.oxygen = np.ones((width, height)) * 0.8       # pO₂ normalizado
        self.glucose = np.ones((width, height)) * 0.7       # Glucosa normalizada
        self.ecm_density = np.zeros((width, height))        # Colágeno + fibronectina
        self.ecm_stiffness = np.zeros((width, height))      # Rigidez (mecanotransducción → YAP)

        # ═══ CITOQUINAS ═══
        self.tgfb = np.zeros((width, height))               # TGF-β
        self.il6 = np.zeros((width, height))                # IL-6 (inflamatoria)
        self.tnfa = np.zeros((width, height))               # TNF-α (inflamatoria)
        self.cxcl12 = np.zeros((width, height))             # CXCL12 / SDF-1 (quimiotaxis)
        self.ifng = np.zeros((width, height))               # IFN-γ (de T-cells activas)

        # ═══ METABOLITOS ═══
        self.lactate = np.zeros((width, height))            # Lactato (Warburg)
        self.ph = np.ones((width, height)) * 7.4            # pH extracelular
        self.ros = np.zeros((width, height))                # ROS extracelulares

        # ═══ FACTORES DE CRECIMIENTO ═══
        self.vegf = np.zeros((width, height))               # VEGF (angiogénesis)
        self.egf = np.ones((width, height)) * 0.1           # EGF basal

        # ═══ VASCULATURA SIMPLIFICADA ═══
        # Vasos en los bordes + algunos vasos aleatorios interiores
        self.vasculature = np.zeros((width, height))
        self.vasculature[0, :] = 1.0
        self.vasculature[-1, :] = 1.0
        self.vasculature[:, 0] = 1.0
        self.vasculature[:, -1] = 1.0
        # Vasos interiores (capilares)
        rng = np.random.default_rng(42)
        n_vessels = max(width * height // 200, 3)
        for _ in range(n_vessels):
            vx, vy = rng.integers(5, width - 5), rng.integers(5, height - 5)
            self.vasculature[max(0, vx-1):min(width, vx+2),
                             max(0, vy-1):min(height, vy+2)] = 0.6

        # ═══ PARÁMETROS: ESCALA FÍSICA Y FICK'S LAW ═══
        self.dx = 20.0  # Escala micrométrica: 1 celda = 20 µm
        
        # Difusión de Oxígeno:
        # Fick's 2nd Law (steady state): C(x) = C0 - (V / 2D) * x^2
        # Con uniform_filter, el operador laplaciano efectivo tiene D equivalente a 1/9.
        # Para que la necrosis (O2=0) ocurra bioquímicamente a ~150 µm de distancia de un capilar:
        # 150 µm = 7.5 celdas. Despejamos V: V = 2 * C0 * D / (x^2) = 2 * 1.0 * (1/9) / (7.5^2) ≈ 0.004.
        self.oxy_supply = 0.025
        self.oxy_consumption_cancer = 0.004  # Calibrado para necrosis a ~150 µm (7.5 celdas)
        self.oxy_consumption_immune = 0.001

        self.glucose_supply = 0.02
        self.glucose_consumption = 0.005
        self.tgfb_decay = 0.05
        self.il6_decay = 0.07
        self.tnfa_decay = 0.08
        self.lactate_decay = 0.04
        self.ecm_decay = 0.008
        self.vegf_decay = 0.06
        self.cxcl12_decay = 0.05
        self.ifng_decay = 0.10

    def update(self, cell_positions: dict):
        """Actualiza todos los campos un paso temporal."""
        cp = cell_positions

        # ═══ OXÍGENO ═══
        self.oxygen = uniform_filter(self.oxygen, size=3, mode='reflect')
        # Suministro desde vasculatura
        supply_mask = self.vasculature > 0
        self.oxygen[supply_mask] = np.minimum(
            self.oxygen[supply_mask] + self.oxy_supply * self.vasculature[supply_mask], 1.0)
        # Consumo por cáncer (alto, Warburg parcial)
        for x, y in cp.get('cancer', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.oxygen[x, y] = max(self.oxygen[x, y] - self.oxy_consumption_cancer, 0.0)
        # Consumo por estroma
        for x, y in cp.get('macrophage', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.oxygen[x, y] = max(self.oxygen[x, y] - self.oxy_consumption_immune, 0.0)

        # ═══ GLUCOSA ═══
        self.glucose = uniform_filter(self.glucose, size=3, mode='reflect')
        supply_mask = self.vasculature > 0
        self.glucose[supply_mask] = np.minimum(
            self.glucose[supply_mask] + self.glucose_supply * self.vasculature[supply_mask], 1.0)
        # Consumo alto por cáncer (Warburg = mucha glucosa)
        for x, y in cp.get('cancer', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.glucose[x, y] = max(self.glucose[x, y] - self.glucose_consumption, 0.0)

        # ═══ LACTATO (producto del efecto Warburg) ═══
        self.lactate = uniform_filter(self.lactate, size=3, mode='reflect')
        for x, y in cp.get('cancer', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                # Más lactato si menos O₂ (Warburg)
                warburg_factor = max(0, 1.0 - self.oxygen[x, y])
                self.lactate[x, y] = min(self.lactate[x, y] + 0.02 * (1 + warburg_factor), 1.0)
        self.lactate *= (1.0 - self.lactate_decay)

        # ═══ pH (baja con lactato) ═══
        self.ph = 7.4 - self.lactate * 1.0  # pH ~6.4-7.4

        # ═══ ROS (aumenta con hipoxia y metabolismo) ═══
        self.ros = np.clip(0.05 + (1.0 - self.oxygen) * 0.3 + self.lactate * 0.15, 0.0, 1.0)

        # ═══ TGF-β ═══
        self.tgfb = uniform_filter(self.tgfb, size=3, mode='reflect')
        for x, y in cp.get('cancer', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.tgfb[x, y] = min(self.tgfb[x, y] + 0.02, 1.0)
        for x, y in cp.get('caf', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.tgfb[x, y] = min(self.tgfb[x, y] + 0.035, 1.0)
        self.tgfb *= (1.0 - self.tgfb_decay)

        # ═══ ECM (secretado por CAFs) ═══
        for x, y in cp.get('caf', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.ecm_density[x, y] = min(self.ecm_density[x, y] + 0.05, 1.0)
                self.ecm_stiffness[x, y] = min(self.ecm_stiffness[x, y] + 0.03, 1.0)
        self.ecm_density *= (1.0 - self.ecm_decay)
        self.ecm_stiffness *= (1.0 - self.ecm_decay * 0.5)  # Rigidez decae más lento
        # ECM → reduce difusión de fármacos (barrera estromal)

        # ═══ IL-6 (proinflamatoria) ═══
        self.il6 = uniform_filter(self.il6, size=3, mode='reflect')
        for x, y in cp.get('caf', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.il6[x, y] = min(self.il6[x, y] + 0.025, 1.0)
        for x, y in cp.get('macrophage', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.il6[x, y] = min(self.il6[x, y] + 0.02, 1.0)
        self.il6 *= (1.0 - self.il6_decay)

        # ═══ TNF-α (inflamatoria, de macrófagos) ═══
        self.tnfa = uniform_filter(self.tnfa, size=3, mode='reflect')
        for x, y in cp.get('macrophage', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.tnfa[x, y] = min(self.tnfa[x, y] + 0.02, 1.0)
        self.tnfa *= (1.0 - self.tnfa_decay)

        # ═══ CXCL12 (quimioquina de CAFs) ═══
        self.cxcl12 = uniform_filter(self.cxcl12, size=3, mode='reflect')
        for x, y in cp.get('caf', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.cxcl12[x, y] = min(self.cxcl12[x, y] + 0.02, 1.0)
        self.cxcl12 *= (1.0 - self.cxcl12_decay)

        # ═══ IFN-γ (de T-cells activas → antitumoral) ═══
        self.ifng = uniform_filter(self.ifng, size=3, mode='reflect')
        for x, y in cp.get('tcell', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                self.ifng[x, y] = min(self.ifng[x, y] + 0.03, 1.0)
        self.ifng *= (1.0 - self.ifng_decay)

        # ═══ VEGF (de células hipóxicas) ═══
        self.vegf = uniform_filter(self.vegf, size=3, mode='reflect')
        for x, y in cp.get('cancer', []):
            if 0 <= x < self.width and 0 <= y < self.height:
                if self.oxygen[x, y] < 0.3:  # Solo en hipoxia
                    self.vegf[x, y] = min(self.vegf[x, y] + 0.03, 1.0)
        self.vegf *= (1.0 - self.vegf_decay)
        # VEGF → neoangiogénesis (crear nuevos vasos)
        high_vegf = self.vegf > 0.6
        self.vasculature[high_vegf] = np.minimum(
            self.vasculature[high_vegf] + 0.005, 0.4)  # Vasos neo, frágiles

        # ═══ CLIP ═══
        self.oxygen = np.clip(self.oxygen, 0.0, 1.0)
        self.glucose = np.clip(self.glucose, 0.0, 1.0)
        self.tgfb = np.clip(self.tgfb, 0.0, 1.0)
        self.ecm_density = np.clip(self.ecm_density, 0.0, 1.0)

    def get_local_conditions(self, x: int, y: int) -> dict:
        """Devuelve todas las condiciones locales del microambiente."""
        if 0 <= x < self.width and 0 <= y < self.height:
            return {
                'oxygen': float(self.oxygen[x, y]),
                'glucose': float(self.glucose[x, y]),
                'tgfb': float(self.tgfb[x, y]),
                'ecm': float(self.ecm_density[x, y]),
                'ecm_stiff': float(self.ecm_stiffness[x, y]),
                'il6': float(self.il6[x, y]),
                'tnfa': float(self.tnfa[x, y]),
                'cxcl12': float(self.cxcl12[x, y]),
                'ifng': float(self.ifng[x, y]),
                'lactate': float(self.lactate[x, y]),
                'ph': float(self.ph[x, y]),
                'ros': float(self.ros[x, y]),
                'vegf': float(self.vegf[x, y]),
                'egf': float(self.egf[x, y]),
            }
        return {k: 0.0 for k in ['oxygen','glucose','tgfb','ecm','ecm_stiff',
                'il6','tnfa','cxcl12','ifng','lactate','ros','vegf','egf','ph']}

    def is_hypoxic(self, x: int, y: int) -> bool:
        if 0 <= x < self.width and 0 <= y < self.height:
            return self.oxygen[x, y] < 0.3
        return False

    def ecm_blocks_movement(self, x: int, y: int, threshold: float = 0.7) -> bool:
        if 0 <= x < self.width and 0 <= y < self.height:
            return self.ecm_density[x, y] > threshold
        return False

    def drug_penetration_factor(self, x: int, y: int) -> float:
        """Factor [0-1] de penetración de fármaco. ECM densa reduce penetración."""
        if 0 <= x < self.width and 0 <= y < self.height:
            ecm = self.ecm_density[x, y]
            return max(1.0 - ecm * 0.6, 0.2)  # Mín 20% penetración
        return 1.0
