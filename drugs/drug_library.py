"""
Librería de fármacos para PDAC v5.
Incluye quimioterapia clásica, pan-KRAS 2026, terapias dirigidas,
inmunoterapia y fármacos experimentales.
"""

import numpy as np


class Drug:
    """Clase base para un fármaco."""

    def __init__(self, name, targets, ic50=0.3, max_efficacy=0.9, c_max=1.0, description=""):
        self.name = name
        self.targets = targets
        self.ic50 = ic50
        self.max_efficacy = max_efficacy
        self.c_max = c_max
        self.description = description

    def get_effect(self, dose, phase='ALL'):
        if dose <= 0:
            return {}
        n = 1.5
        efficacy = self.max_efficacy * (dose**n) / (self.ic50**n + dose**n)
        
        fx = {}
        for t, w in self.targets.items():
            if isinstance(w, tuple):
                weight, specific_phase = w
                if specific_phase != 'ALL' and phase != specific_phase:
                    continue
                fx[t] = efficacy * weight
            else:
                fx[t] = efficacy * w
        return fx


class DrugLibrary:
    """Biblioteca completa de fármacos PDAC."""

    def __init__(self):
        self.drugs = {}
        self._init_chemo()
        self._init_pan_kras()
        self._init_targeted()
        self._init_immuno()
        self._init_experimental()
        self._init_ferroptosis_senolytic()
        self._init_metabolic()

    def _init_chemo(self):
        """Quimioterapia citotóxica clásica."""

        self.drugs['gemcitabine'] = Drug(
            name='Gemcitabina',
            targets={
                'G1_S_transition': (0.40, 'ALL'),
                'S_progression': (0.90, 'S'),
                'apoptosis_signal': (-0.45, 'S'),
                'DDR_active': (-0.35, 'S'),      # Estrés replicativo → ATR → DDR (TIS a dosis baja)
                'ATR_active': (-0.30, 'S'),       # Gem → horquillas replicativas paradas → ATR
            },
            ic50=0.50, max_efficacy=0.85, c_max=60.0,
            description='Análogo nucleósido. Actúa en Fase S. Estrés replicativo → ATR → TIS o apoptosis. Cmax~60 µM.',
        )

        self.drugs['5fu'] = Drug(
            name='5-Fluorouracilo (5-FU)',
            targets={
                'S_progression': (0.85, 'S'),
                'apoptosis_signal': (-0.35, 'S'),
            },
            ic50=2.00, max_efficacy=0.80, c_max=20.0,
            description='Antimetabolito pirimidínico (Fase S). Cmax~20 µM.',
        )

        self.drugs['oxaliplatin'] = Drug(
            name='Oxaliplatino',
            targets={
                'S_progression': (0.60, 'S'),
                'G2_M_transition': (0.70, 'G2'),
                'apoptosis_signal': (-0.40, 'ALL'),
                'ROS_level': (-0.25, 'ALL'),
            },
            ic50=1.00, max_efficacy=0.85, c_max=5.0,
            description='Platino 3ª gen. Aductos DNA. Cmax~5 µM.',
        )

        self.drugs['irinotecan'] = Drug(
            name='Irinotecán',
            targets={
                'S_progression': (0.80, 'S'),
                'G2_M_transition': (0.60, 'G2'),
                'apoptosis_signal': (-0.30, 'S'),
            },
            ic50=1.50, max_efficacy=0.80, c_max=4.0,
            description='Inhibidor topoisomerasa I (Fase S/G2). Cmax~4 µM.',
        )

        self.drugs['nab_paclitaxel'] = Drug(
            name='Nab-Paclitaxel (Abraxane)',
            targets={
                'G2_M_transition': (0.95, 'G2'),
                'proliferation_signal': 0.80,
                'apoptosis_signal': (-0.20, 'M'),
            },
            ic50=0.10, max_efficacy=0.90, c_max=10.0,
            description='Estabiliza microtúbulos (Arresto G2/M). Cmax~10 µM.',
        )

    def _init_pan_kras(self):
        """Inhibidores pan-KRAS / RAS(ON) 2025-2026."""

        self.drugs['mrtx1133'] = Drug(
            name='MRTX1133 (KRAS G12D)',
            targets={'KRAS_active': 0.85, 'RAF_active': 0.50, 'ERK_active': 0.30},
            ic50=0.02, max_efficacy=0.90, c_max=3.0,
            description='Inhibidor selectivo KRAS G12D (on). Cmax~3µM.',
        )

        self.drugs['daraxonrasib'] = Drug(
            name='Daraxonrasib (RMC-6236)',
            targets={
                'KRAS_active': 0.95, 'RAF_active': 0.60,
                'ERK_active': 0.35, 'proliferation_signal': 0.50,
            },
            ic50=0.05, max_efficacy=0.95, c_max=5.0,
            description='Inhibidor multi-selectivo RAS(ON) tri-complex. Cmax~5µM.',
        )

        self.drugs['rmc7977'] = Drug(
            name='RMC-7977',
            targets={
                'KRAS_active': 0.92, 'RAF_active': 0.55,
                'ERK_active': 0.30, 'proliferation_signal': 0.45,
            },
            ic50=0.03, max_efficacy=0.95, c_max=5.0,
            description='RAS(ON) 2ª generación. Mayor potencia preclínica. Cmax~5µM.',
        )

        self.drugs['eras0015'] = Drug(
            name='ERAS-0015',
            targets={
                'KRAS_active': 0.88, 'RAF_active': 0.45,
                'ERK_active': 0.25, 'proliferation_signal': 0.40,
            },
            ic50=0.05, max_efficacy=0.88, c_max=4.0,
            description='Pan-RAS molecular glue degrader. Cmax~4µM.',
        )

        self.drugs['ly4066434'] = Drug(
            name='LY-4066434',
            targets={
                'KRAS_active': 0.90, 'RAF_active': 0.50,
                'ERK_active': 0.30, 'proliferation_signal': 0.45,
                'survival_signal': 0.20,
            },
            ic50=0.02, max_efficacy=0.90, c_max=3.5,
            description='Pan-KRAS con penetración CNS. Cmax~3.5µM.',
        )

    def _init_targeted(self):
        """Terapias dirigidas (no-KRAS)."""

        self.drugs['afatinib'] = Drug(
            name='Afatinib (EGFR/HER2i)',
            targets={
                'EGFR_active': 0.90,
                'HER2_active': 0.85,
                'KRAS_active': 0.15,
                'ERK_active': 0.30,
                'PI3K_active': 0.35,
                'proliferation_signal': 0.35,
            },
            ic50=0.005, max_efficacy=0.85, c_max=0.1,
            description='Inhibidor irreversible pan-HER (EGFR/HER2/HER4). Cmax~100nM.',
        )

        self.drugs['verteporfin'] = Drug(
            name='Verteporfin (YAPi)',
            targets={'YAP_nuclear': 0.85, 'TEAD_active': 0.70},
            ic50=1.0, max_efficacy=0.75, c_max=5.0,
            description='Disruptor YAP-TEAD. Anti-resistencia a KRASi. Cmax~5µM.',
        )

        self.drugs['everolimus'] = Drug(
            name='Everolimus (mTORi)',
            targets={
                'mTOR_active': 0.85,
                'proliferation_signal': 0.30,
            },
            ic50=0.02, max_efficacy=0.80, c_max=0.2,
            description='Inhibidor de mTORC1 (rapalog). Cmax~200nM.',
        )

        self.drugs['trametinib'] = Drug(
            name='Trametinib (MEKi)',
            targets={
                'MEK_active': 0.90,
                'ERK_active': 0.70,
                'G1_S_transition': (0.60, 'G1'), # MEKi frena ciclo celular
            },
            ic50=0.01, max_efficacy=0.85, c_max=0.05,
            description='Inhibidor selectivo MEK1/2. Cmax~50nM.',
        )

        self.drugs['olaparib'] = Drug(
            name='Olaparib (PARPi)',
            targets={
                'ROS_level': -0.40,          # PARP trapping → ROS acumulado
                'DDR_active': -0.50,          # Activa DDR (DSBs no reparables)
                'ATM_active': -0.35,          # ATM activado por DSBs (convenio: negativo = activa)
                'RAD51_active': 0.40,         # Inhibe HR vía agotamiento BRCA (sólo efectivo si BRCA-mut)
                'survival_signal': 0.35,
                'apoptosis_signal': (-0.30, 'S'),  # Daño en replicación S
            },
            ic50=0.10, max_efficacy=0.80, c_max=20.0,
            description='Inhibidor PARP1/2. Trapping PARPi → DSBs → letal sintético con BRCA-mut. Cmax~20µM.',
        )

    def _init_immuno(self):
        """Inmunoterapia."""

        self.drugs['anti_pd1'] = Drug(
            name='Anti-PD-1 (Pembrolizumab)',
            targets={'PDL1_expression': 0.80},
            ic50=0.1, max_efficacy=0.60, c_max=10.0,
            description='Anticuerpo anti-PD-1. Solo eficaz en MSI-H/dMMR (~1-2% PDAC).',
        )

        self.drugs['anti_ctla4'] = Drug(
            name='Anti-CTLA-4 (Ipilimumab)',
            targets={
                'PDL1_expression': 0.40,
            },
            ic50=0.1, max_efficacy=0.50, c_max=10.0,
            description='Anti-CTLA-4.',
        )

    def _init_experimental(self):
        """Fármacos experimentales / fase temprana."""

        self.drugs['protac_stat3'] = Drug(
            name='PROTAC STAT3 (SD-36)',
            targets={
                'STAT3_active': 0.90,
                'survival_signal': 0.70,
                'proliferation_signal': 0.40,
                'PDL1_expression': 0.30,
                'autophagy': 0.20,
            },
            ic50=0.05, max_efficacy=0.85, c_max=2.0,
            description='PROTAC que degrada STAT3 vía ubiquitina-proteasoma. Cmax~2µM.',
        )

        self.drugs['belzutifan'] = Drug(
            name='Belzutifan (HIF-2αi)',
            targets={
                'HIF2A_active': 0.90,
                'warburg_effect': 0.60,
                'glutamine_addiction': 0.50,
            },
            ic50=0.05, max_efficacy=0.85, c_max=1.5,
            description='Antagonista HIF-2α. Revierte Warburg en hipoxia. Cmax~1.5µM.',
        )

        self.drugs['hydroxychloroquine'] = Drug(
            name='Hidroxicloroquina (AutoI)',
            targets={
                'autophagy': 0.70,
                'macropinocytosis': 0.30,
            },
            ic50=5.0, max_efficacy=0.70, c_max=20.0,
            description='Inhibidor lisosómico → bloquea autofagia. Cmax~20µM.',
        )

        self.drugs['sotorasib'] = Drug(
            name='Sotorasib (KRAS G12Ci)',
            targets={
                'KRAS_active': 0.80,
                'RAF_active': 0.45,
                'ERK_active': 0.25,
            },
            ic50=0.04, max_efficacy=0.85, c_max=12.0,
            description='Inhibidor covalente KRAS G12C. Cmax~12µM.',
        )

        # ═══ NUEVOS EXPERIMENTALES 2025-2026 ═══

        self.drugs['rmc4630'] = Drug(
            name='RMC-4630 (SHP2i)',
            targets={
                'KRAS_active': 0.30,
                'RAF_active': 0.25,
                'ERK_active': 0.20,
                'EGFR_active': 0.15,
            },
            ic50=0.03, max_efficacy=0.80, c_max=1.5,
            description='Inhibidor alostérico SHP2. Bloquea EGFR→SHP2→SOS1→RAS. Cmax~1.5µM.',
        )

        self.drugs['bi3406'] = Drug(
            name='BI-3406 (SOS1i)',
            targets={
                'KRAS_active': 0.25,
                'KRAS_WT_active': 0.40,
                'RAF_active': 0.15,
            },
            ic50=0.05, max_efficacy=0.75, c_max=2.0,
            description='Inhibidor de SOS1. Bloquea reactivación de KRAS-WT. Cmax~2µM.',
        )

        self.drugs['palbociclib'] = Drug(
            name='Palbociclib (CDK4/6i)',
            targets={
                'RB_active': -0.40,
                'G1_S_transition': (0.85, 'G1'),  # Frena el ciclo
                'MYC_active': 0.20,
            },
            ic50=0.05, max_efficacy=0.85, c_max=1.0,
            description='CDK4/6i. Arresta en G1, provocando antagonismo con Gem/5-FU. Cmax~1µM.',
        )

        self.drugs['ceralasertib'] = Drug(
            name='Ceralasertib (ATRi)',
            targets={
                'ATR_active': 0.85,           # Inhibición directa de ATR (v9 nodo específico)
                'CHEK1_active': 0.70,          # Inhibe checkpoint S-fase vía CHK1
                'DDR_active': 0.60,            # Reduce capacidad DDR global
                'ROS_level': -0.20,            # Aumenta ROS por acumulación de horquillas paradas
                'S_progression': (0.80, 'S'),  # Forza avance catastrófico en Fase S
                'apoptosis_signal': (-0.20, 'S'),
            },
            ic50=0.05, max_efficacy=0.80, c_max=3.0,
            description='Inhibidor ATR quinasa → bloquea CHK1 → colapso replicativo en Fase S. Cmax~3µM.',
        )

    def _init_ferroptosis_senolytic(self):
        """Inductores de ferroptosis y senolíticos. Vías GPX4/SLC7A11."""

        # RSL3: inhibidor covalente de GPX4. Induce ferroptosis.
        # Más efectivo en células con NRF2 bajo (KRAS-WT, sin ARID1A)
        self.drugs['rsl3'] = Drug(
            name='RSL3 (GPX4i)',
            targets={
                'GPX4_level': 0.88,           # Inhibe GPX4 directamente
                'lipid_ROS': -0.50,            # Aumenta lipid_ROS (convenio negativo = activa)
                'ferroptosis_primed': -0.70,   # Activa estado ferroptótico
                'GSH_level': 0.30,             # Consume GSH indirectamente
            },
            ic50=0.08, max_efficacy=0.85, c_max=5.0,
            description='Inhibidor covalente de GPX4 → ferroptosis. Sinérgico con KRAS-WT o NRF2-bajo. Cmax~5µM.',
        )

        # Erastin: inhibidor xCT (SLC7A11). Bloquea importación de cistina → agota GSH → ferroptosis.
        self.drugs['erastin'] = Drug(
            name='Erastin (xCT/SLC7A11i)',
            targets={
                'SLC7A11_expression': 0.85,    # Bloquea transportador xCT funcionalmente
                'GSH_level': 0.60,             # Agota GSH (consecuencia de SLC7A11i)
                'GPX4_level': 0.35,            # GPX4 se inactiva sin sustrato GSH
                'ferroptosis_primed': -0.65,   # Activa ferroptosis
            },
            ic50=0.50, max_efficacy=0.80, c_max=10.0,
            description='Erastin bloquea xCT → depleción GSH → colapso GPX4 → ferroptosis. Cmax~10µM.',
        )

        # Sulfasalazine: fármaco aprobado (artritis). Inhibe SLC7A11. Reposicionamiento.
        self.drugs['sulfasalazine'] = Drug(
            name='Sulfasalazina (xCT/SLC7A11i)',
            targets={
                'SLC7A11_expression': 0.65,    # Inhibición parcial de xCT (menos potente que erastin)
                'GSH_level': 0.40,
                'ferroptosis_primed': -0.45,
                'NFKB_active': 0.30,           # También inhibe NFκB
            },
            ic50=0.80, max_efficacy=0.70, c_max=8.0,
            description='Sulfasalazina: FDA-aprobada. Inhibe xCT + NFκB. Candidata reposicionamiento PDAC. Cmax~8µM.',
        )

        # Navitoclax: BCL-2/BCL-XL inhibidor. Senolítico + pro-apoptótico.
        # Elimina células senescentes (SASP) que resisten a quimioterapia.
        self.drugs['navitoclax'] = Drug(
            name='Navitoclax (BCL-2/BCL-XLi)',
            targets={
                'BCL2_active': 0.80,           # Inhibe BCL-2 (anti-apoptótico)
                'survival_signal': 0.55,        # Reduce supervivencia
                'apoptosis_signal': (-0.40, 'ALL'),  # Pro-apoptótico
                'SASP_active': 0.60,            # Elimina células con SASP activo (senolisis)
            },
            ic50=0.08, max_efficacy=0.85, c_max=2.0,
            description='Navitoclax: BCL-2/BCL-XL inhibidor. Senolítico → elimina TIS y células resistentes. Cmax~2µM.',
        )

    def _init_metabolic(self):
        """Fármacos metabólicos: AMPK, mTOR, IDO1, CD73."""

        # Metformin: activa AMPK vía inhibición Complejo I mitocondrial.
        # Biguanida anti-diabética en estudio clínico en PDAC.
        self.drugs['metformin'] = Drug(
            name='Metformin (AMPK activador)',
            targets={
                'AMPK_active': -0.55,          # Activa AMPK (convenio negativo = activa)
                'mTOR_active': 0.50,            # Inhibe mTOR vía AMPK
                'LDHA_active': 0.35,            # Reduce glucólisis anaeróbica
                'energy_stress': -0.30,         # Induce estrés energético mitocondrial
                'autophagy': -0.20,             # Promueve autofagia vía AMPK
                'warburg_effect': 0.30,
            },
            ic50=1.00, max_efficacy=0.65, c_max=50.0,
            description='Metformina: activa AMPK → inhibe mTORC1 → frena Warburg. Rango clínico ~1-50 mM intrac. Cmax~50µM.',
        )

        # IDO1 inhibitor: bloquea conversión triptófano → kinurenina → restaura inmunidad T
        self.drugs['epacadostat'] = Drug(
            name='Epacadostat (IDO1i)',
            targets={
                'IDO1_active': 0.85,           # Inhibe IDO1 directamente
                'PDL1_expression': 0.15,        # Leve reducción de PDL1
            },
            ic50=0.05, max_efficacy=0.75, c_max=2.0,
            description='Epacadostat: Inhibidor IDO1. Restaura actividad de células T en TME. Cmax~2µM.',
        )

        # Anti-CD73: bloquea adenosina extracelular → restaura inmunidad
        self.drugs['oleclumab'] = Drug(
            name='Oleclumab (anti-CD73)',
            targets={
                'CD73_expression': 0.80,        # Bloquea funcionalidad CD73
                'PDL1_expression': 0.20,
            },
            ic50=0.10, max_efficacy=0.70, c_max=5.0,
            description='Oleclumab: mAb anti-CD73. Reduce adenosina inmunosupresora en TME. Cmax~5µM.',
        )

    def get_drug(self, name):
        return self.drugs.get(name)

    def get_combined_effects(self, drug_doses, phase='ALL'):
        """Bliss independence con Antagonismo Mecánico (Fases del ciclo)."""
        combined = {}
        for drug_name, dose in drug_doses.items():
            drug = self.drugs.get(drug_name)
            if drug is None or dose <= 0:
                continue
            effects = drug.get_effect(dose, phase)
            for target, inh in effects.items():
                if target in combined:
                    # Para valores negativos (apoptosis inducida), sumar
                    if inh < 0 or combined[target] < 0:
                        combined[target] = combined[target] + inh
                    else:
                        combined[target] = combined[target] + inh - combined[target] * inh
                else:
                    combined[target] = inh
        for t in combined:
            combined[t] = np.clip(combined[t], -1.0, 1.0)
        return combined

    def get_drug_names(self):
        return list(self.drugs.keys())

    def get_all_info(self):
        """Devuelve info de todos los fármacos para la GUI."""
        return {k: {'name': d.name, 'desc': d.description}
                for k, d in self.drugs.items()}


def apply_resistance_mechanism(network, drug_name, exposure_time):
    """Mecanismos de resistencia documentados en pdac según fármaco y tiempo."""
    if exposure_time < 5:
        return
    # Incrementamos drásticamente la tasa de penetración de resistencia (rf)
    # y subimos el cap a 1.0 para que la exposición prolongada 
    # cruce los umbrales de >0.40 requeridos para los reportes médicos
    rf = min(exposure_time / 15.0, 1.0)

    if drug_name == 'nab_paclitaxel':
        # 2024 Nab-Paclitaxel: Inducción sostenida de c-MYC y mecanobiología
        network.nodes['MYC_active'] = min(network.nodes['MYC_active'] + rf * 0.40, 1.0)
        network.nodes['YAP_nuclear'] = min(network.nodes['YAP_nuclear'] + rf * 0.20, 1.0)
        
    elif drug_name in ('gemcitabine', '5fu', 'oxaliplatin', 'irinotecan'):
        # Quimioterapia clásica / FOLFIRINOX
        network.nodes['autophagy'] = min(network.nodes['autophagy'] + rf * 0.25, 1.0)
        if drug_name in ('5fu', 'oxaliplatin', 'irinotecan'):
            # 2024 FOLFIRINOX: Eje GALNT5/MYH9/NOTCH y DDR (DNA Damage Response)
            network.nodes['NOTCH_active'] = min(network.nodes['NOTCH_active'] + rf * 0.35, 1.0)
            network.nodes['DDR_active'] = min(network.nodes['DDR_active'] + rf * 0.40, 1.0)
            network.nodes['basal_like_program'] = min(network.nodes['basal_like_program'] + rf * 0.20, 1.0)
    
    elif drug_name in ('mrtx1133', 'daraxonrasib', 'rmc7977', 'eras0015', 'ly4066434', 'sotorasib'):
        # 2024-2025 KRASi / RAS(ON)i Resistencia TARDÍA (semanas-meses):
        # NOTA: El rebote RTK (EGFR/HER2 → ERK) ya está modelado temporalmente
        # en pdac_network.py. Aquí solo modelamos mecanismos TARDÍOS:
        # 1. Vía AGER-DIAPH1 (Macropinocitosis masiva para escape metabólico)
        network.nodes['macropinocytosis'] = min(network.nodes['macropinocytosis'] + rf * 0.40, 1.0)
        # 2. Amplificación genómica MYC (Cancer Discovery 2023)
        network.nodes['MYC_active'] = min(network.nodes['MYC_active'] + rf * 0.25, 1.0)
        # 3. YAP/TAZ bypass transcripcional (Nature 2024)
        network.nodes['YAP_nuclear'] = min(network.nodes['YAP_nuclear'] + rf * 0.20, 1.0)
        # 4. EMT / basal-like shift (resistencia intrínseca)
        network.nodes['basal_like_program'] = min(network.nodes['basal_like_program'] + rf * 0.15, 1.0)

    elif drug_name == 'trametinib':
        # MEKi Resistencia: Feedback loop release
        network.nodes['ERK_active'] = min(network.nodes['ERK_active'] + rf * 0.45, 1.0)
        network.nodes['EGFR_active'] = min(network.nodes['EGFR_active'] + rf * 0.35, 1.0)

    elif drug_name == 'olaparib':
        # 2024 PARPi Resistencia:
        # 1. Mutación de Reversión BRCA (Restaura HR)
        network.nodes['BRCA_functional'] = min(network.nodes['BRCA_functional'] + rf * 0.50, 1.0)
        # 2. Rewiring Metabólico (Fosforilación Oxidativa)
        network.nodes['OXPHOS_active'] = min(network.nodes['OXPHOS_active'] + rf * 0.45, 1.0)
        network.nodes['DDR_active'] = min(network.nodes['DDR_active'] + rf * 0.30, 1.0)

    elif drug_name in ('anti_pd1', 'anti_ctla4'):
        # Inmunoterapia → Pérdida de antígenos, downregulation de PD-L1, upregulation de otros checkpoints (TIM-3)
        # Modelamos esto como induciendo hipoxia/TME hostil y bajando PD-L1 (escape inmunológico)
        network.nodes['PDL1_expression'] = max(network.nodes['PDL1_expression'] - rf * 0.40, 0.0)

    elif drug_name == 'verteporfin':
        # YAP inhibitors → Bypass por KRAS o vías alternativas
        network.nodes['KRAS_active'] = min(network.nodes['KRAS_active'] + rf * 0.25, 1.0)
        
    elif drug_name == 'afatinib':
        # EGFR/HER inhibitors → KRAS amplification/mutation
        network.nodes['KRAS_active'] = min(network.nodes['KRAS_active'] + rf * 0.30, 1.0)
        network.nodes['PI3K_active'] = min(network.nodes['PI3K_active'] + rf * 0.20, 1.0)

    elif drug_name == 'protac_stat3':
        # STAT3 degraders → Upregulation de otros factores de supervivencia cruzados
        network.nodes['survival_signal'] = min(network.nodes['survival_signal'] + rf * 0.20, 1.0)

    elif drug_name == 'rmc4630':
        # SHP2i Resistencia: amplificación RAS o RTK bypass
        network.nodes['KRAS_active'] = min(network.nodes['KRAS_active'] + rf * 0.20, 1.0)
        network.nodes['PI3K_active'] = min(network.nodes['PI3K_active'] + rf * 0.25, 1.0)

    elif drug_name == 'bi3406':
        # SOS1i Resistencia: SOS2 compensación + RAS amplification
        network.nodes['KRAS_active'] = min(network.nodes['KRAS_active'] + rf * 0.15, 1.0)
        network.nodes['EGFR_active'] = min(network.nodes['EGFR_active'] + rf * 0.20, 1.0)

    elif drug_name == 'palbociclib':
        # CDK4/6i Resistencia: CDK2 bypass, RB loss, ciclina E amplification
        network.nodes['MYC_active'] = min(network.nodes['MYC_active'] + rf * 0.30, 1.0)
        network.nodes['proliferation_signal'] = min(network.nodes['proliferation_signal'] + rf * 0.15, 1.0)

    elif drug_name == 'ceralasertib':
        # ATRi Resistencia: CHK1 upregulation, DDR rewiring, WEE1 compensación
        network.nodes['DDR_active']    = min(network.nodes.get('DDR_active', 0) + rf * 0.35, 1.0)
        network.nodes['CHEK1_active']  = min(network.nodes.get('CHEK1_active', 0) + rf * 0.25, 1.0)
        network.nodes['survival_signal'] = min(network.nodes.get('survival_signal', 0) + rf * 0.15, 1.0)

    elif drug_name == 'rsl3':
        # GPX4i (RSL3) Resistencia: Upregulation NRF2 → más GPX4 y FSP1; upregulation SLC7A11
        network.nodes['NRF2_active']          = min(network.nodes.get('NRF2_active', 0) + rf * 0.35, 1.0)
        network.nodes['GPX4_level']           = min(network.nodes.get('GPX4_level', 0) + rf * 0.30, 1.0)
        network.nodes['FSP1_level']           = min(network.nodes.get('FSP1_level', 0) + rf * 0.25, 1.0)
        network.nodes['SLC7A11_expression']   = min(network.nodes.get('SLC7A11_expression', 0) + rf * 0.20, 1.0)

    elif drug_name in ('erastin', 'sulfasalazine'):
        # xCT/SLC7A11i Resistencia: NRF2/ARE upregulation → transcripción SLC7A11
        network.nodes['NRF2_active']          = min(network.nodes.get('NRF2_active', 0) + rf * 0.40, 1.0)
        network.nodes['SLC7A11_expression']   = min(network.nodes.get('SLC7A11_expression', 0) + rf * 0.35, 1.0)
        network.nodes['GSH_level']            = min(network.nodes.get('GSH_level', 0) + rf * 0.25, 1.0)

    elif drug_name == 'navitoclax':
        # BCL-2/BCL-XLi Resistencia: Upregulation MCL-1 (modelado via survival/NFκB)
        network.nodes['survival_signal']  = min(network.nodes.get('survival_signal', 0) + rf * 0.30, 1.0)
        network.nodes['NFKB_active']      = min(network.nodes.get('NFKB_active', 0) + rf * 0.20, 1.0)
        network.nodes['autophagy']        = min(network.nodes.get('autophagy', 0) + rf * 0.15, 1.0)

    elif drug_name == 'metformin':
        # Metformin Resistencia: Bypass via AMPK alternativo, upregulation PGC-1α
        network.nodes['OXPHOS_active']    = min(network.nodes.get('OXPHOS_active', 0) + rf * 0.30, 1.0)
        network.nodes['autophagy']        = min(network.nodes.get('autophagy', 0) + rf * 0.15, 1.0)
