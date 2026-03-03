"""
Modelo principal del tumor PDAC v7.
- PK temporal (Cmax, t½)
- Tregs, MDSCs, NK cells
- CAF subtypes tracking
- Enhanced history
Compatible con Mesa 3.x API.
"""

import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import mesa
from agents.cancer_cell import CancerCell
from agents.caf import CAF
from agents.immune_cells import MacrophageM2, CD8_Tcell, Treg, MDSC, NKCell
from microenvironment.diffusion import Microenvironment
from drugs.drug_library import DrugLibrary, apply_resistance_mechanism


# ═══ FARMACOCINÉTICA SIMPLIFICADA ═══
DRUG_PK = {
    # {drug: (t_half_hours, t_max_hours, dosing_interval_hours)}
    # dosing_interval = período entre dosis (simula infusiones repetidas)
    'gemcitabine': (1.0, 0.5, 8),     # IV, t½ ~1h, re-dosis cada 8h (sim)
    '5fu': (0.5, 0.3, 6),             # IV continua, re-dosis cada 6h (sim)
    'oxaliplatin': (16.0, 2.0, 48),   # IV q2w, re-dosis cada 48h
    'irinotecan': (6.0, 1.5, 24),     # IV, re-dosis cada 24h
    'nab_paclitaxel': (20.0, 1.0, 48),# IV semanal
    'daraxonrasib': (6.0, 2.0, 12),   # Oral BID
    'rmc7977': (8.0, 3.0, 24),
    'eras0015': (5.0, 2.0, 12),
    'ly4066434': (7.0, 2.5, 24),
    'mrtx1133': (4.0, 1.5, 8),
    'sotorasib': (5.0, 1.0, 24),      # Oral QD
    'afatinib': (37.0, 3.0, 24),      # Oral QD
    'verteporfin': (5.0, 2.0, 12),
    'trametinib': (120.0, 1.5, 24),   # Oral QD
    'everolimus': (30.0, 1.0, 24),    # Oral QD
    'olaparib': (12.0, 1.5, 12),      # Oral BID
    'anti_pd1': (600.0, 24.0, 336),   # IV q2w
    'anti_ctla4': (360.0, 24.0, 504), # IV q3w
    'protac_stat3': (8.0, 2.0, 12),
    'belzutifan': (14.0, 2.0, 24),    # Oral QD
    'hydroxychloroquine': (40.0, 3.0, 24),
}


# ═══ FARMACOCINÉTICA INTRACELULAR ═══
# Metabolitos activos intracelulares: persisten MUCHO más que el fármaco plasmático.
# Esta es la causa raíz del fallo de validación clínica: gemcitabine plasma t½=1h,
# pero dFdCTP intracelular t½=25h → actividad sostenida entre dosis semanales.
DRUG_INTRACELL_PK = {
    # {drug: (t_half_intracell_hours, uptake_rate, phosphorylation_eff)}
    # t_half_intracell: vida media del metabolito activo intracelular
    # uptake_rate: tasa de incorporación plasmático → intracelular [0-1]
    # phosphorylation_eff: eficiencia de fosforilación a metabolito activo [0-1]
    'gemcitabine':    (25.0, 0.40, 0.70),  # dFdCTP t½=25h (Heinemann 2003)
    '5fu':            (2.0,  0.35, 0.60),  # FUTP/FdUMP t½=2h
    'oxaliplatin':    (72.0, 0.30, 0.85),  # Aductos ADN t½=72h (Ehrsson 2002)
    'irinotecan':     (4.0,  0.25, 0.55),  # SN-38 activo t½=4h
    'olaparib':       (8.0,  0.45, 0.90),  # PARPi intracelular t½=8h
    'nab_paclitaxel': (12.0, 0.35, 0.80),  # Paclitaxel intracelular t½=12h
}


def pk_concentration(t_hours, t_half, t_max, dose, dosing_interval=None):
    """
    Modelo PK 1-compartimento con dosificación repetida.
    Si dosing_interval > 0, simula dosis periódicas (sawtooth).
    Returns: concentración en [0, dose].
    """
    if t_hours <= 0:
        return 0.0

    # Dosificación periódica: usar tiempo dentro del ciclo actual
    if dosing_interval and dosing_interval > 0:
        t_hours = (t_hours % dosing_interval) or dosing_interval

    if t_hours < t_max:
        return dose * (t_hours / t_max)
    else:
        ke = np.log(2) / t_half
        return dose * np.exp(-ke * (t_hours - t_max))


class TumorModel(mesa.Model):
    """Gemelo digital PDAC v8 — multi-escala clínico."""

    def __init__(self, width=50, height=50, n_cancer=30, n_caf=15,
                 n_macrophage=10, n_tcell=8, n_treg=0, n_mdsc=0, n_nk=0,
                 seed=None, mutation_profile=None):
        super().__init__(rng=seed)
        self.width = width
        self.height = height
        
        self.current_step = 0
        self._rng = np.random.default_rng(seed)
        
        if mutation_profile is None:
            from signaling.pdac_network import generate_mutation_profile
            self.mutation_profile = generate_mutation_profile(self._rng)
        else:
            self.mutation_profile = mutation_profile
        
        self.msi_status = self.mutation_profile.get('MSI_status', 'MSS')

        self.grid = mesa.space.SingleGrid(width, height, torus=False)
        self.microenv = Microenvironment(width, height)

        self.drug_library = DrugLibrary()
        self.drug_doses = {}
        self.drug_effects = {}
        self.active_drugs = []
        self.pk_factor = 1.0  # Factor PK [0-1]

        # Concentraciones intracelulares acumuladas (metabolitos activos)
        # Clave: drug_name → concentración normalizada [0-1] en células tumorales
        self.intracell_conc = {}   # se actualiza cada step con cinética 1-compartimento

        # Tregs/MDSCs auto-escalados si no especificados
        if n_treg == 0:
            n_treg = max(n_tcell // 3, 2)
        if n_mdsc == 0:
            n_mdsc = max(n_macrophage // 3, 2)
        if n_nk == 0:
            n_nk = max(n_tcell // 4, 1)

        # ═══ CARRYING CAPACITY (Gompertz) ═══
        self.carrying_capacity = int(width * height * 0.75)
        self.growth_inhibition = 0.0  # [0-1] factor de inhibición por densidad

        # ═══ PROTOCOLO CLÍNICO ═══
        self.protocol_name = None
        self._protocol_active = False

        # ═══ RECLUTAMIENTO INMUNE DINÁMICO ═══
        self._base_tcell_rate = n_tcell  # Para calcular tasa de reclutamiento
        self._base_mac_rate = n_macrophage

        # ═══ NECROTIC CORE ═══
        self._hypoxia_tracker = {}  # {agent_id: horas en hipoxia severa}

        # Historial enriquecido
        self.history = {
            'step': [], 'cancer_alive': [], 'cancer_dead': [],
            'resistant_count': [], 'resistant_pct': [],
            'caf_count': [], 'mycaf_count': [], 'icaf_count': [], 'apcaf_count': [],
            'macrophage_count': [], 'mac_m1_pct': [],
            'tcell_count': [], 'tcell_exhaustion': [],
            'treg_count': [], 'mdsc_count': [], 'nk_count': [],
            'avg_oxygen': [], 'avg_lactate': [], 'avg_ph': [],
            'avg_proliferation': [], 'avg_apoptosis': [],
            'basal_like_pct': [], 'avg_pdl1': [],
            'pk_factor': [], 'ca199': [],
            'necrotic_count': [], 'growth_inhibition': [],
            # v9: nuevas métricas
            'senescent_count': [],        # Células en senescencia activa
            'avg_ferroptosis_risk': [],   # Riesgo medio de ferroptosis
            'avg_nrf2': [],               # Actividad media NRF2 (defensa antioxidante)
            'avg_ampk': [],               # Actividad media AMPK
        }
        self.initial_node_averages = {}  # Para Reportes Delta T0 vs Tf

        self._create_agents(n_cancer, n_caf, n_macrophage, n_tcell,
                            n_treg, n_mdsc, n_nk)

    def _create_agents(self, n_cancer, n_caf, n_macrophage, n_tcell,
                        n_treg, n_mdsc, n_nk):
        cx, cy = self.width // 2, self.height // 2

        for _ in range(n_cancer):
            pos = self._find_empty_near(cx, cy, 8)
            if pos is None: break
            a = CancerCell(self, _pos=pos, mutations=self.mutation_profile)
            self.grid.place_agent(a, pos)

        for _ in range(n_caf):
            pos = self._find_empty_near(cx, cy, 14)
            if pos is None: break
            a = CAF(self, pos=pos)
            a.activated = True
            a.activation_level = 0.5
            self.grid.place_agent(a, pos)

        for _ in range(n_macrophage):
            pos = self._find_empty_random()
            if pos is None: break
            self.grid.place_agent(MacrophageM2(self, pos=pos), pos)

        for _ in range(n_tcell):
            pos = self._find_empty_border()
            if pos is None: break
            self.grid.place_agent(CD8_Tcell(self, pos=pos), pos)

        for _ in range(n_treg):
            pos = self._find_empty_random()
            if pos is None: break
            self.grid.place_agent(Treg(self, pos=pos), pos)

        for _ in range(n_mdsc):
            pos = self._find_empty_random()
            if pos is None: break
            self.grid.place_agent(MDSC(self, pos=pos), pos)

        for _ in range(n_nk):
            pos = self._find_empty_border()
            if pos is None: break
            self.grid.place_agent(NKCell(self, pos=pos), pos)

    def _find_empty_near(self, cx, cy, radius):
        for _ in range(200):
            x = int(np.clip(cx + self._rng.integers(-radius, radius+1), 0, self.width-1))
            y = int(np.clip(cy + self._rng.integers(-radius, radius+1), 0, self.height-1))
            if self.grid.is_cell_empty((x, y)):
                return (x, y)
        return None

    def _find_empty_random(self):
        for _ in range(200):
            x, y = int(self._rng.integers(0, self.width)), int(self._rng.integers(0, self.height))
            if self.grid.is_cell_empty((x, y)):
                return (x, y)
        return None

    def _find_empty_border(self):
        for _ in range(200):
            side = int(self._rng.integers(4))
            if side == 0:   pos = (0, int(self._rng.integers(0, self.height)))
            elif side == 1: pos = (self.width-1, int(self._rng.integers(0, self.height)))
            elif side == 2: pos = (int(self._rng.integers(0, self.width)), 0)
            else:           pos = (int(self._rng.integers(0, self.width)), self.height-1)
            if self.grid.is_cell_empty(pos):
                return pos
        return None

    def set_protocol(self, protocol_name):
        """Activar protocolo clínico (scheduling automático)."""
        self.protocol_name = protocol_name
        self._protocol_active = True

    def set_drug_doses(self, doses):
        self.drug_doses = doses
        self.active_drugs = [n for n, d in doses.items() if d > 0]
        tmp_uM = {}
        for dn in self.active_drugs:
            d = doses[dn]
            drug = self.drug_library.get_drug(dn)
            c_max = drug.c_max if drug else 1.0
            tmp_uM[dn] = d * c_max
        self.current_doses_uM = tmp_uM
        
        # Efectos pre-calculados por fase
        self.drug_effects_by_phase = {
            ph: self.drug_library.get_combined_effects(tmp_uM, phase=ph) 
            for ph in ('G1', 'S', 'G2', 'M', 'ALL')
        }
        self.drug_effects = self.drug_effects_by_phase['ALL']  # Legacy fallback

    def step(self):
        self.current_step += 1
        
        # ═══ DELTA TRACKING (T0) ═══
        # Capturamos el estado ESTABILIZADO de la red (tras warm-up)
        # para que los deltas porcentuales reflejen el cambio real vs baseline.
        if self.current_step == 1:
            # Warm-up: correr la red 5 ciclos sin fármacos para estabilizar
            cancer_cells = [a for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
            for cell in cancer_cells:
                if hasattr(cell, 'network'):
                    for _ in range(5):
                        cell.network.update(microenv={}, drug_effects={}, active_drugs=[])
            self.initial_node_averages = self._get_average_signaling()

        # ═══ PROTOCOLO CLÍNICO: scheduling automático ═══
        if self._protocol_active and self.protocol_name:
            from drugs.treatment_protocols import get_active_drugs_for_hour
            doses = get_active_drugs_for_hour(self.protocol_name, self.current_step)
            if doses:
                self.drug_doses = doses
                self.active_drugs = [n for n, d in doses.items() if d > 0]
            else:
                self.drug_doses = {}
                self.active_drugs = []

        # ═══ PK: concentración plasmática en µM ═══
        plasma_uM = {}
        if self.active_drugs:
            for drug_name in self.active_drugs:
                dose_fraction = self.drug_doses.get(drug_name, 0)
                if dose_fraction <= 0:
                    continue
                drug = self.drug_library.get_drug(drug_name)
                c_max = drug.c_max if drug else 1.0
                pk_data = DRUG_PK.get(drug_name, (12.0, 1.0, 24))
                t_half, t_max, d_interval = pk_data
                conc_fraction = pk_concentration(self.current_step, t_half, t_max, 1.0, d_interval)
                plasma_uM[drug_name] = conc_fraction * c_max * dose_fraction
            self.pk_factor = 1.0  # Legacy support
        else:
            self.pk_factor = 1.0

        # ═══ PK INTRACELULAR: metabolitos activos ═══
        # Para fármacos con metabolitos activos intracelulares (gemcitabine→dFdCTP, etc.),
        # la concentración efectiva en la célula depende tanto de la concentración plasmática
        # actual como de la acumulación previa (t½ intracelular mucho mayor).
        # Modelo: dC_ic/dt = uptake * C_plasma - (ln2/t½_ic) * C_ic
        # Discretizado: C_ic(t+1) = C_ic(t) * exp(-ke_ic * Δt) + uptake * C_plasma(t) * Δt
        dt = 1.0  # 1 hora por step
        for drug_name, pk_ic in DRUG_INTRACELL_PK.items():
            t_half_ic, uptake_rate, phospho_eff = pk_ic
            ke_ic = np.log(2) / t_half_ic
            c_plasma = plasma_uM.get(drug_name, 0.0)
            c_ic_prev = self.intracell_conc.get(drug_name, 0.0)
            # Cinética: decaimiento intracelular + captación desde plasma
            c_ic_new = c_ic_prev * np.exp(-ke_ic * dt) + uptake_rate * phospho_eff * c_plasma * dt
            self.intracell_conc[drug_name] = max(c_ic_new, 0.0)

        # Concentración efectiva = max(plasma, intracelular) para fármacos con PK_IC
        # Para el resto, usar plasma directamente
        self.current_doses_uM = {}
        for drug_name in self.active_drugs:
            if drug_name in DRUG_INTRACELL_PK and self.intracell_conc.get(drug_name, 0) > 0:
                # Usar el mayor de plasma e intracelular (intracelular domina cuando acumula)
                self.current_doses_uM[drug_name] = max(
                    plasma_uM.get(drug_name, 0.0),
                    self.intracell_conc[drug_name]
                )
            else:
                self.current_doses_uM[drug_name] = plasma_uM.get(drug_name, 0.0)
        # También incluir fármacos intracelulares que persisten aunque no estén activos ahora
        for drug_name, c_ic in self.intracell_conc.items():
            if drug_name not in self.current_doses_uM and c_ic > 0.01:
                self.current_doses_uM[drug_name] = c_ic

        # Pre-calcular efectos combinados por fase del ciclo (Ahorra CPU y permite Antagonismo real)
        self.drug_effects_by_phase = {
             'G1':  self.drug_library.get_combined_effects(self.current_doses_uM, phase='G1'),
             'S':   self.drug_library.get_combined_effects(self.current_doses_uM, phase='S'),
             'G2':  self.drug_library.get_combined_effects(self.current_doses_uM, phase='G2'),
             'M':   self.drug_library.get_combined_effects(self.current_doses_uM, phase='M'),
             'ALL': self.drug_library.get_combined_effects(self.current_doses_uM, phase='ALL')
        }
        self.drug_effects = self.drug_effects_by_phase['ALL']

        # ═══ CARRYING CAPACITY (Gompertz) ═══
        alive_cancer = sum(1 for a in self.agents
                          if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive)
        density = alive_cancer / max(self.carrying_capacity, 1)
        self.growth_inhibition = min(density ** 1.5, 0.95)  # Gompertz-like

        # Posiciones por tipo
        cell_pos = {'cancer': [], 'caf': [], 'macrophage': [], 'tcell': []}
        for a in self.agents:
            if not hasattr(a, 'cell_type') or a.pos is None:
                continue
            ct = a.cell_type
            if ct == 'cancer':     cell_pos['cancer'].append(a.pos)
            elif ct == 'caf':      cell_pos['caf'].append(a.pos)
            elif ct == 'macrophage': cell_pos['macrophage'].append(a.pos)
            elif ct in ('tcell', 'nk'): cell_pos['tcell'].append(a.pos)

        # Microambiente
        self.microenv.update(cell_pos)

        # ═══ NECROTIC CORE (hipoxia prolongada) ═══
        for a in list(self.agents):
            if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive and a.pos:
                ox = self.microenv.oxygen[a.pos[0], a.pos[1]]
                aid = a.unique_id
                if ox < 0.05:
                    self._hypoxia_tracker[aid] = self._hypoxia_tracker.get(aid, 0) + 1
                    if self._hypoxia_tracker[aid] > 24:  # >24h en hipoxia severa
                        a.alive = False
                        if a.pos is not None:
                            self.grid.remove_agent(a)
                else:
                    self._hypoxia_tracker[aid] = max(self._hypoxia_tracker.get(aid, 0) - 1, 0)

        # Resistencia (drug library)
        for a in list(self.agents):
            if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive and hasattr(a, 'signaling'):
                for dn, exp in a.drug_exposure.items():
                    apply_resistance_mechanism(a.signaling, dn, exp)

        # Step todos los agentes
        self.agents.shuffle_do("step")

        # ═══ RECLUTAMIENTO INMUNE DINÁMICO ═══
        if self.current_step % 24 == 0:  # Cada 24h
            self._recruit_immune_cells()

        self._collect()

    def _collect(self):
        cancer = [a for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer']
        alive = [c for c in cancer if c.alive]
        n = max(len(alive), 1)

        avg_p = avg_a = avg_pdl1 = avg_ferrop = avg_nrf2 = avg_ampk = 0.0
        basal = resistant = senescent_n = 0
        for c in alive:
            if hasattr(c, 'signaling'):
                avg_p += c.signaling.get_proliferation_rate()
                avg_a += c.signaling.get_apoptosis_probability()
                avg_pdl1 += c.signaling.get_pdl1_level()
                avg_ferrop += c.signaling.get_ferroptosis_risk()
                avg_nrf2   += c.signaling.nodes.get('NRF2_active', 0.0)
                avg_ampk   += c.signaling.nodes.get('AMPK_active', 0.0)
                if c.signaling.is_basal_like():
                    basal += 1
            if getattr(c, 'resistant', False):
                resistant += 1
            if getattr(c, 'senescent', False):
                senescent_n += 1

        # CAF subtypes
        cafs = [a for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'caf']
        mycaf = sum(1 for c in cafs if getattr(c, 'subtype', '') == 'myCAF')
        icaf = sum(1 for c in cafs if getattr(c, 'subtype', '') == 'iCAF')
        apcaf = sum(1 for c in cafs if getattr(c, 'subtype', '') == 'apCAF')

        # Macrophage polarization
        macs = [a for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'macrophage']
        m1_count = sum(1 for m in macs if getattr(m, 'polarization', 'M2') == 'M1')
        m1_pct = (m1_count / max(len(macs), 1)) * 100

        # T-cell exhaustion
        tcells = [a for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'tcell' and a.alive]
        avg_exhaust = np.mean([t.exhaustion for t in tcells]) if tcells else 0.0

        # Others
        tregs = sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'treg')
        mdscs = sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'mdsc')
        nks = sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'nk')

        h = self.history
        h['step'].append(self.current_step)
        h['cancer_alive'].append(len(alive))
        h['cancer_dead'].append(len(cancer) - len(alive))
        h['resistant_count'].append(resistant)
        h['resistant_pct'].append(resistant / n * 100)
        h['caf_count'].append(len(cafs))
        h['mycaf_count'].append(mycaf)
        h['icaf_count'].append(icaf)
        h['apcaf_count'].append(apcaf)
        h['macrophage_count'].append(len(macs))
        h['mac_m1_pct'].append(m1_pct)
        h['tcell_count'].append(len(tcells))
        h['tcell_exhaustion'].append(float(avg_exhaust))
        h['treg_count'].append(tregs)
        h['mdsc_count'].append(mdscs)
        h['nk_count'].append(nks)
        h['avg_oxygen'].append(float(np.mean(self.microenv.oxygen)))
        h['avg_lactate'].append(float(np.mean(self.microenv.lactate)))
        h['avg_ph'].append(float(np.mean(self.microenv.ph)))
        h['avg_proliferation'].append(avg_p / n)
        h['avg_apoptosis'].append(avg_a / n)
        h['basal_like_pct'].append(basal / n * 100)
        h['avg_pdl1'].append(avg_pdl1 / n)
        h['pk_factor'].append(float(self.pk_factor))
        h['senescent_count'].append(senescent_n)
        h['avg_ferroptosis_risk'].append(avg_ferrop / n)
        h['avg_nrf2'].append(avg_nrf2 / n)
        h['avg_ampk'].append(avg_ampk / n)

        # CA 19-9 (biomarcador sérico)
        from simulation.clinical_endpoints import estimate_ca199
        baseline = h['cancer_alive'][0] if h['cancer_alive'] else 30
        ca199 = estimate_ca199(len(alive), baseline, self.mutation_profile)
        h['ca199'].append(float(ca199))

        # Necrotic count
        necrotic = sum(1 for a in self.agents
                       if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and not a.alive)
        h['necrotic_count'].append(necrotic)
        h['growth_inhibition'].append(float(self.growth_inhibition))

    def _recruit_immune_cells(self):
        """
        Reclutamiento inmune dinámico: cada 24h, nuevas células inmunes
        llegan desde la vasculatura. Tasa proporcional a señales inflamatorias.
        """
        # Señales que atraen inmunidad
        avg_il6 = float(np.mean(self.microenv.il6))
        avg_ifng = float(np.mean(self.microenv.ifng))
        inflammation = (avg_il6 + avg_ifng) / 2.0

        # Reclutamiento base + inflamación
        n_new_tc = max(1, int(self._base_tcell_rate * 0.1 * (1 + inflammation * 2)))
        n_new_mac = max(1, int(self._base_mac_rate * 0.08 * (1 + inflammation)))

        # Limitar para no saturar
        n_new_tc = min(n_new_tc, 3)
        n_new_mac = min(n_new_mac, 2)

        for _ in range(n_new_tc):
            pos = self._find_empty_border()
            if pos:
                self.grid.place_agent(CD8_Tcell(self, pos=pos), pos)

        for _ in range(n_new_mac):
            pos = self._find_empty_random()
            if pos:
                self.grid.place_agent(MacrophageM2(self, pos=pos), pos)

    def get_grid_state(self):
        state = np.zeros((self.width, self.height), dtype=int)
        for a in self.agents:
            if a.pos is None or not hasattr(a, 'cell_type'):
                continue
            x, y = a.pos
            ct = a.cell_type
            if ct == 'cancer':
                state[x, y] = 1 if a.alive else 2
            elif ct == 'caf':        state[x, y] = 3
            elif ct == 'macrophage': state[x, y] = 4
            elif ct == 'tcell':      state[x, y] = 5
            elif ct == 'treg':       state[x, y] = 6
            elif ct == 'mdsc':       state[x, y] = 7
            elif ct == 'nk':         state[x, y] = 8
        return state

    def get_history_dataframe(self):
        return pd.DataFrame(self.history)

    def get_summary(self):
        alive = sum(1 for a in self.agents
                    if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive)
        return {
            'Paso': self.current_step,
            'Células tumorales': alive,
            'CAFs': sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'caf'),
            'Macrófagos M2': sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'macrophage'),
            'T-cells CD8+': sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'tcell' and a.alive),
            'Tregs': sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'treg'),
            'MDSCs': sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'mdsc'),
            'NK cells': sum(1 for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'nk'),
            'O₂ promedio': f"{np.mean(self.microenv.oxygen):.2f}",
        }

    def _get_average_signaling(self):
        """Devuelve el promedio de actividad de todas las proteínas en el tumor vivo."""
        cancer_cells = [a for a in self.agents if hasattr(a, 'cell_type') and a.cell_type == 'cancer' and a.alive]
        if not cancer_cells:
            return {}
            
        avg = {}
        for c in cancer_cells:
            if hasattr(c, 'signaling'):
                for k, v in c.signaling.nodes.items():
                    avg[k] = avg.get(k, 0) + v
        n = len(cancer_cells)
        for k in avg:
            avg[k] /= n
        return avg
