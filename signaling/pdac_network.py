"""
Red de señalización molecular PDAC v9 — Modelo molecular expandido.
~67 nodos | 15 vías | Ferroptosis + Senescencia + AMPK/LKB1 + IDO1/CD73 + DDR detallado.

Nuevas vías añadidas (v9):
  - Ferroptosis: GPX4, SLC7A11/xCT, GSH, lipid_ROS, FSP1, NRF2 (KEAP1-NRF2 axis)
  - AMPK/LKB1: sensor energético, antagonista mTORC1, activa autofagia
  - Metabolismo expandido: LDHA, GLS, FASN (lipogénesis)
  - DDR detallado: ATM, ATR, CHEK1, CHEK2, RAD51 (HR vs NHEJ)
  - Inmunoevasión ampliada: IDO1 (triptófano), CD73 (adenosina), CD47 (fagocitosis)
  - Senescencia: SASP_active (para célula senescente)

Compatibilidad total con v5.1: todos los nodos antiguos se conservan.
Basado en literatura PDAC 2024-2025: TCGA, ICGC, NCI PDAC Atlas.
"""

import numpy as np


# ═══════════════════════════════════════════════════════════
#  PERFILES KRAS (efectos biológicos diferenciados por variante)
# ═══════════════════════════════════════════════════════════
KRAS_VARIANTS = {
    'G12D': {
        'kras_level': 0.90,
        'pi3k_boost': 0.35,
        'erk_boost': 0.20,
        'macro_boost': 0.40,
        'nrf2_boost': 0.20,
        'description': 'Más frecuente (~41%). Fuerte PI3K/AKT, ferroptosis-resistente vía NRF2.',
    },
    'G12V': {
        'kras_level': 0.95,
        'pi3k_boost': 0.20,
        'erk_boost': 0.35,
        'macro_boost': 0.30,
        'nrf2_boost': 0.25,
        'description': 'Agresivo (~30%). Alta ERK, proliferación intensa, NRF2 elevado.',
    },
    'G12R': {
        'kras_level': 0.85,
        'pi3k_boost': 0.15,
        'erk_boost': 0.15,
        'macro_boost': 0.70,
        'nrf2_boost': 0.15,
        'description': 'Intermedio (~16%). Muy dependiente macropinocitosis y glutaminolysis.',
    },
    'G12C': {
        'kras_level': 0.80,
        'pi3k_boost': 0.15,
        'erk_boost': 0.20,
        'macro_boost': 0.30,
        'nrf2_boost': 0.15,
        'description': 'Raro en PDAC (~1-2%). Sensible a Sotorasib.',
    },
    'Q61H': {
        'kras_level': 0.92,
        'pi3k_boost': 0.25,
        'erk_boost': 0.30,
        'macro_boost': 0.35,
        'nrf2_boost': 0.20,
        'description': 'Raro (~2%). Resistente a RAS(ON) inhibitors.',
    },
    'WT': {
        'kras_level': 0.10,
        'pi3k_boost': 0.05,
        'erk_boost': 0.05,
        'macro_boost': 0.10,
        'nrf2_boost': 0.05,
        'description': 'Wild-type (~8%). EGFR-dependiente. NRF2 bajo → más sensible a ferroptosis.',
    },
}


class NodeDict(dict):
    """
    Diccionario interceptor que aplica efectos farmacológicos on-the-fly.
    Convención: inh > 0 → inhibe (value *= 1-inh) | inh < 0 → activa (value -= inh)
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.inh_map = {}

    def __setitem__(self, key, value):
        if key in self.inh_map:
            inh = self.inh_map[key]
            if inh < 0:
                value = min(value - inh, 1.0)
            else:
                value = max(value * (1.0 - inh), 0.0)
        super().__setitem__(key, value)


class SignalingNetwork:
    """
    Red de señalización molecular PDAC v9.
    ~67 nodos, 15 vías, calibrada para dinámica realista (DT ~24-72h).
    """

    def __init__(self, mutations: dict = None, seed: int = None):
        self.mutations = mutations or {}
        self.kras_variant = self.mutations.get('KRAS_variant', 'G12D')
        self.kras_profile = KRAS_VARIANTS.get(self.kras_variant, KRAS_VARIANTS['G12D'])
        self.drug_exposure_steps = 0

        rng = np.random.default_rng(seed)
        self.noise_prolif = 1.0 + rng.normal(0, 0.12)
        self.noise_surv   = 1.0 + rng.normal(0, 0.08)
        self.noise_apop   = 1.0 + rng.normal(0, 0.10)

        self.nodes = NodeDict({
            # RAS-MAPK
            'EGFR_active': 0.25, 'HER2_active': 0.10,
            'KRAS_active': 0.10, 'KRAS_WT_active': 0.05,
            'NRAS_active': 0.05, 'RAF_active': 0.10,
            'MEK_active': 0.10,  'ERK_active': 0.10,
            # PI3K/AKT/mTOR
            'PI3K_active': 0.10, 'PTEN_active': 0.80,
            'AKT_active': 0.10,  'TSC_active': 0.70,
            'mTOR_active': 0.10, 'S6K_active': 0.10,
            # AMPK/LKB1 (nuevo v9)
            'LKB1_active': 0.80, 'AMPK_active': 0.20, 'energy_stress': 0.10,
            # YAP/TAZ
            'LATS_active': 0.60, 'YAP_nuclear': 0.10,
            'TEAD_active': 0.10, 'CTGF_expression': 0.05, 'CYR61_expression': 0.05,
            # Hipoxia
            'oxygen_level': 0.80, 'HIF1A_active': 0.00,
            'HIF2A_active': 0.00, 'VEGF_expression': 0.05,
            # TGF-β/SMAD
            'TGFB_signal': 0.10, 'SMAD4_active': 0.10, 'SMAD23_active': 0.10,
            # WNT
            'WNT_active': 0.10, 'BCAT_nuclear': 0.10,
            # NOTCH
            'NOTCH_active': 0.20, 'HES1_expression': 0.10,
            # NFκB
            'NFKB_active': 0.15,
            # JAK/STAT
            'IL6_signal': 0.10, 'STAT3_active': 0.15,
            # MYC
            'MYC_active': 0.20,
            # Supresores
            'TP53_active': 0.80, 'P21_active': 0.50,
            'CDKN2A_active': 0.80, 'RB_active': 0.70,
            # Apoptosis
            'BCL2_active': 0.50, 'BAX_active': 0.20, 'CASP3_active': 0.05,
            # Metabolismo expandido (v9)
            'autophagy': 0.20,        'macropinocytosis': 0.10,
            'warburg_effect': 0.10,   'glutamine_addiction': 0.10,
            'lipogenesis': 0.10,      'OXPHOS_active': 0.20,
            'LDHA_active': 0.30,      'GLS_active': 0.20,     'FASN_active': 0.20,
            # NRF2/KEAP1 (nuevo v9)
            'NRF2_active': 0.20,
            # Ferroptosis (nuevo v9)
            'SLC7A11_expression': 0.30, 'GSH_level': 0.55,
            'GPX4_level': 0.60,         'FSP1_level': 0.25,
            'lipid_ROS': 0.10,          'ferroptosis_primed': 0.00,
            # EMT
            'SNAIL_active': 0.05, 'ZEB1_active': 0.05,
            'ECAD_expression': 0.80, 'VIM_expression': 0.10,
            'basal_like_program': 0.00,
            # DDR detallado (v9)
            'DDR_active': 0.30,        'ATM_active': 0.10,    'ATR_active': 0.10,
            'CHEK1_active': 0.10,      'CHEK2_active': 0.10,
            'BRCA_functional': 0.90,   'RAD51_active': 0.40,  'ROS_level': 0.15,
            # Inmunoevasión ampliada (v9)
            'PDL1_expression': 0.10,   'MHC1_expression': 0.70,
            'CXCL12_secretion': 0.10,  'IDO1_active': 0.05,
            'CD73_expression': 0.15,   'CD47_expression': 0.25,
            # Senescencia (v9)
            'SASP_active': 0.00,
            # Salidas (compatibilidad v5.1)
            'proliferation_signal': 0.30, 'G1_S_transition': 0.10,
            'S_progression': 0.10,        'G2_M_transition': 0.10,
            'survival_signal': 0.50,      'apoptosis_signal': 0.10,
            'invasion_signal': 0.05,      'angiogenesis_signal': 0.05,
        })

        self._apply_mutations()

    def _apply_mutations(self):
        m = self.mutations

        if m.get('KRAS') == 'HET_MUT':
            kp = self.kras_profile
            self.nodes['KRAS_active']     = kp['kras_level']
            self.nodes['macropinocytosis'] = kp['macro_boost']
            self.nodes['autophagy']        = max(0.3, kp['macro_boost'] * 0.6)
            self.nodes['NRF2_active']      = min(0.20 + kp['nrf2_boost'], 0.80)
        else:
            self.nodes['KRAS_active'] = 0.10

        tp53 = m.get('TP53', 'WT')
        if tp53 == 'HOM_LOSS':
            self.nodes['TP53_active'] = 0.01
            self.nodes['P21_active']  = 0.05
            self.nodes['BAX_active']  = 0.05
            self.nodes['DDR_active']  = 0.20
            self.nodes['SLC7A11_expression'] = min(self.nodes['SLC7A11_expression'] + 0.20, 1.0)
        elif tp53 == 'HET_LOSS':
            self.nodes['TP53_active'] = 0.40
            self.nodes['P21_active']  = 0.30
            self.nodes['BAX_active']  = 0.15
            self.nodes['DDR_active']  = 0.25

        cdkn2a = m.get('CDKN2A', 'WT')
        if cdkn2a == 'HOM_LOSS':
            self.nodes['CDKN2A_active'] = 0.01
            self.nodes['RB_active']     = 0.20
        elif cdkn2a == 'HET_LOSS':
            self.nodes['CDKN2A_active'] = 0.40
            self.nodes['RB_active']     = 0.50

        smad4 = m.get('SMAD4', 'WT')
        if smad4 == 'HOM_LOSS':
            self.nodes['SMAD4_active']  = 0.01
            self.nodes['SMAD23_active'] = 0.05
        elif smad4 == 'HET_LOSS':
            self.nodes['SMAD4_active']  = 0.40
            self.nodes['SMAD23_active'] = 0.30

        if m.get('YAP_amp') == 'AMP':
            self.nodes['YAP_nuclear'] = 0.75
            self.nodes['LATS_active'] = 0.20

        brca = m.get('BRCA', 'WT')
        if brca == 'HOM_LOSS':
            self.nodes['BRCA_functional'] = 0.05
            self.nodes['RAD51_active']    = 0.10
            self.nodes['DDR_active']      = 0.10
        elif brca == 'HET_LOSS':
            self.nodes['BRCA_functional'] = 0.40
            self.nodes['RAD51_active']    = 0.25
            self.nodes['DDR_active']      = 0.20

        pten = m.get('PTEN', 'WT')
        if pten == 'HOM_LOSS':
            self.nodes['PTEN_active'] = 0.02
        elif pten == 'HET_LOSS':
            self.nodes['PTEN_active'] = 0.40

        if m.get('MYC') == 'AMP':
            self.nodes['MYC_active'] = 0.80

        lkb1 = m.get('LKB1', 'WT')
        if lkb1 == 'HOM_LOSS':
            self.nodes['LKB1_active'] = 0.02
        elif lkb1 == 'HET_LOSS':
            self.nodes['LKB1_active'] = 0.40
        else:
            self.nodes['LKB1_active'] = 0.80

        if m.get('ARID1A', 'WT') in ('HOM_LOSS', 'HET_LOSS'):
            self.nodes['NRF2_active'] = min(self.nodes['NRF2_active'] + 0.25, 1.0)

    def update(self, microenv: dict = None, drug_effects: dict = None,
               active_drugs: list = None, cell_cycle_phase: str = 'G1'):
        me = microenv or {}
        de = drug_effects or {}
        active_drugs = active_drugs or []
        kp = self.kras_profile
        variant = self.kras_variant

        # Especificidad farmacológica KRAS
        inh_map = {}
        kras_specific_drugs = {
            'mrtx1133':    ['G12D'],
            'sotorasib':   ['G12C'],
            'daraxonrasib':['G12D', 'G12V', 'G12C'],
            'rmc7977':     ['G12D', 'G12V', 'G12C', 'G12R'],
            'eras0015':    ['G12D', 'G12V', 'G12C', 'G12R', 'Q61H'],
            'ly4066434':   ['G12D', 'G12V', 'G12C', 'G12R', 'Q61H'],
        }
        for target, inhibition in de.items():
            actual_inh = inhibition
            if target == 'KRAS_active' and inhibition > 0:
                is_mismatch = False
                for drug_id in active_drugs:
                    allowed = kras_specific_drugs.get(drug_id)
                    if allowed and variant not in allowed:
                        is_mismatch = True
                        break
                if is_mismatch:
                    actual_inh = 0.0 if variant == 'WT' else actual_inh * 0.1
            inh_map[target] = actual_inh

        self.nodes.inh_map = inh_map
        n = self.nodes

        # Tracking temporal KRASi
        kras_drug_names = set(kras_specific_drugs.keys())
        has_krasi = bool(kras_drug_names & set(active_drugs))
        if has_krasi and inh_map.get('KRAS_active', 0) > 0.3:
            self.drug_exposure_steps += 1
        else:
            self.drug_exposure_steps = max(0, self.drug_exposure_steps - 2)
        temporal_factor = min(self.drug_exposure_steps / 48.0, 1.0)

        self._apply_mutations()

        # ─── TME INPUTS ───
        if 'oxygen' in me:
            n['oxygen_level'] = np.clip(me['oxygen'], 0.0, 1.0)
        oxy = n['oxygen_level']
        glucose = me.get('glucose', 0.7)

        if 'tgfb' in me:
            n['TGFB_signal'] = np.clip(n['TGFB_signal'] * 0.7 + me['tgfb'] * 0.3, 0.0, 1.0)
            if me['tgfb'] > 0.3:
                n['VIM_expression']     = min(n['VIM_expression'] + me['tgfb'] * 0.02, 1.0)
                n['ECAD_expression']    = max(n['ECAD_expression'] - me['tgfb'] * 0.01, 0.0)
                n['basal_like_program'] = min(n['basal_like_program'] + me['tgfb'] * 0.01, 1.0)
        if 'il6' in me:
            n['IL6_signal'] = np.clip(n['IL6_signal'] * 0.7 + me['il6'] * 0.3, 0.0, 1.0)
            if me['il6'] > 0.2:
                n['STAT3_active'] = min(n['STAT3_active'] + me['il6'] * 0.03, 1.0)
                n['NFKB_active']  = min(n['NFKB_active']  + me['il6'] * 0.02, 1.0)
        if 'ecm_stiff' in me and me['ecm_stiff'] > 0.4:
            n['YAP_nuclear'] = min(n['YAP_nuclear'] + me['ecm_stiff'] * 0.03, 1.0)
        if 'egf' in me:
            n['EGFR_active'] = np.clip(n['EGFR_active'] * 0.6 + me['egf'] * 0.4, 0.0, 1.0)
        if glucose < 0.3:
            n['autophagy']        = min(n['autophagy'] + 0.05, 1.0)
            n['macropinocytosis'] = min(n['macropinocytosis'] + 0.03, 1.0)
        if me.get('lactate', 0) > 0.3:
            n['PDL1_expression'] = min(n['PDL1_expression'] + 0.02, 1.0)
            n['warburg_effect']  = min(n['warburg_effect'] + me['lactate'] * 0.02, 1.0)
        ifng = me.get('ifng', 0)
        if ifng > 0.2:
            n['MHC1_expression'] = min(n['MHC1_expression'] + ifng * 0.05, 1.0)
            n['PDL1_expression'] = min(n['PDL1_expression'] + ifng * 0.03, 1.0)
            n['IDO1_active']     = min(n['IDO1_active'] + ifng * 0.04, 1.0)

        # ─── RAS-MAPK ───
        if self.mutations.get('KRAS') != 'HET_MUT':
            n['KRAS_active'] = np.clip(n['EGFR_active'] * 0.5 + n['HER2_active'] * 0.5 + 0.15, 0.0, 0.8)
        kras_total = np.clip(n['KRAS_active'] + n['KRAS_WT_active'] * 0.15, 0.0, 1.0)
        raf_in = kras_total * 0.85
        n['RAF_active'] = np.clip(raf_in, 0.0, 1.0) if kras_total > 0.1 else n['RAF_active'] * 0.6
        n['MEK_active'] = np.clip(n['RAF_active'] * 0.80, 0.0, 1.0) if n['RAF_active'] > 0.2 else n['MEK_active'] * 0.6
        erk_in = n['MEK_active'] * 0.80 + kp['erk_boost'] * kras_total * 0.3
        n['ERK_active'] = np.clip(erk_in, 0.0, 1.0) if n['MEK_active'] > 0.2 else n['ERK_active'] * 0.6
        if n['ERK_active'] > 0.6:
            n['RAF_active'] *= (1.0 - (n['ERK_active'] - 0.6) * 0.3)

        # RTK rebound (Barbacid PNAS 2025)
        if n['ERK_active'] < 0.25:
            rebound_maturity = min(max(self.drug_exposure_steps - 6, 0) / 42.0, 1.0)
            rtk_rebound = (0.25 - n['ERK_active']) * 1.2 * rebound_maturity
            egfr_new = min(n['EGFR_active'] + rtk_rebound, 0.85)
            her2_new = min(n['HER2_active'] + rtk_rebound * 0.8, 0.75)
            wt_new   = min(n['KRAS_WT_active'] + rtk_rebound * 0.5, 0.60)
            if 'EGFR_active' in inh_map and inh_map['EGFR_active'] > 0.2:
                n['EGFR_active'] = egfr_new
            else:
                dict.__setitem__(n, 'EGFR_active', egfr_new)
            if 'HER2_active' in inh_map and inh_map['HER2_active'] > 0.2:
                n['HER2_active'] = her2_new
            else:
                dict.__setitem__(n, 'HER2_active', her2_new)
            dict.__setitem__(n, 'KRAS_WT_active', wt_new)

        # ─── AMPK/LKB1 (v9) — se calcula ANTES de mTOR para inhibición continua ───
        energy_stress_val = max(0.0, 1.0 - glucose * 0.65 - oxy * 0.35)
        n['energy_stress'] = np.clip(energy_stress_val, 0.0, 1.0)
        ampk_in = n['LKB1_active'] * 0.35 + n['energy_stress'] * 0.30 + n['autophagy'] * 0.10
        n['AMPK_active'] = np.clip(ampk_in, 0.0, 1.0)
        if n['AMPK_active'] > 0.15:
            n['autophagy'] = min(n['autophagy'] + n['AMPK_active'] * 0.04, 1.0)

        # ─── PI3K/AKT/mTOR (mTOR integra inhibición AMPK continuamente) ───
        pi3k_in = (kras_total * 0.4 + kp['pi3k_boost'] * kras_total * 0.3
                   + n['EGFR_active'] * 0.3 + n['HER2_active'] * 0.3)
        pi3k_in *= (1.0 - n['PTEN_active'] * 0.5)
        n['PI3K_active'] = np.clip(pi3k_in, 0.0, 1.0)
        n['AKT_active']  = np.clip(n['PI3K_active'] * 0.75, 0.0, 1.0) if n['PI3K_active'] > 0.2 else n['AKT_active'] * 0.6
        # AMPK inhibe mTORC1 vía TSC2 y RAPTOR (Shaw 2004) — relación continua, no condicional
        mtor_in = (n['AKT_active'] * 0.7
                   * (1.0 - n['TSC_active'] * 0.3)
                   * (1.0 - n['AMPK_active'] * 0.40))
        n['mTOR_active'] = np.clip(mtor_in, 0.0, 1.0)
        n['S6K_active']  = np.clip(n['mTOR_active'] * 0.6, 0.0, 1.0)
        n['BCL2_active'] = np.clip(0.3 + n['AKT_active'] * 0.4, 0.0, 1.0)

        # ─── YAP/TAZ ───
        if n['LATS_active'] > 0.5:
            n['YAP_nuclear'] *= (1.0 - (n['LATS_active'] - 0.5) * 0.2)
        if n['YAP_nuclear'] > 0.3:
            n['TEAD_active']      = np.clip(n['YAP_nuclear'] * 0.85, 0.0, 1.0)
            n['CTGF_expression']  = np.clip(n['TEAD_active'] * 0.6,  0.0, 1.0)
            n['CYR61_expression'] = np.clip(n['TEAD_active'] * 0.5,  0.0, 1.0)
        else:
            n['TEAD_active'] = n['YAP_nuclear'] * 0.3

        # ─── HIPOXIA/HIF ───
        if oxy < 0.2:
            n['HIF1A_active'] = np.clip(1.0 - oxy * 3, 0.0, 1.0)
            n['HIF2A_active'] = np.clip(0.7 - oxy * 2, 0.0, 1.0)
        elif oxy < 0.4:
            n['HIF1A_active'] = np.clip(0.5 - oxy, 0.0, 1.0)
            n['HIF2A_active'] = np.clip(0.3 - oxy * 0.5, 0.0, 1.0)
        else:
            n['HIF1A_active'] *= 0.4
            n['HIF2A_active'] *= 0.4
        hif = max(n['HIF1A_active'], n['HIF2A_active'])
        if hif > 0.3:
            n['warburg_effect']      = np.clip(hif * 0.7, 0.0, 1.0)
            n['glutamine_addiction'] = np.clip(hif * 0.5, 0.0, 1.0)
            n['VEGF_expression']     = np.clip(hif * 0.6, 0.0, 1.0)
            n['angiogenesis_signal'] = n['VEGF_expression'] * 0.7
            n['PDL1_expression']     = min(n['PDL1_expression'] + hif * 0.03, 1.0)
            n['MHC1_expression']     = max(n['MHC1_expression'] - 0.02, 0.1)
            n['CD73_expression']     = min(n['CD73_expression'] + hif * 0.02, 1.0)
        else:
            n['warburg_effect']      = max(n['warburg_effect'] * 0.8, 0.05)
            n['glutamine_addiction'] = max(n['glutamine_addiction'] * 0.8, 0.03)

        # ─── TGF-β/SMAD ───
        if n['TGFB_signal'] > 0.2:
            n['SMAD23_active'] = np.clip(n['TGFB_signal'] * 0.6, 0.0, 1.0)
        if n['TGFB_signal'] > 0.3 and n['SMAD4_active'] > 0.3:
            n['apoptosis_signal'] = min(n['apoptosis_signal'] + 0.03, 0.20)
        if n['TGFB_signal'] > 0.3 and n['SMAD4_active'] < 0.2:
            n['SNAIL_active']       = min(n['SNAIL_active'] + 0.04, 1.0)
            n['ZEB1_active']        = min(n['ZEB1_active'] + 0.03, 1.0)
            n['basal_like_program'] = min(n['basal_like_program'] + 0.03, 1.0)

        # ─── EMT ───
        emt = (n['SNAIL_active'] + n['ZEB1_active']) / 2.0
        if emt > 0.2:
            n['ECAD_expression'] = max(n['ECAD_expression'] - emt * 0.05, 0.0)
            n['VIM_expression']  = min(n['VIM_expression']  + emt * 0.05, 1.0)
            n['invasion_signal'] = np.clip(emt * 0.5, 0.0, 1.0)

        # ─── WNT ───
        if n['WNT_active'] > 0.3:
            n['BCAT_nuclear'] = np.clip(n['WNT_active'] * 0.6, 0.0, 1.0)

        # ─── NOTCH ───
        if n['NOTCH_active'] > 0.4:
            n['HES1_expression']    = np.clip(n['NOTCH_active'] * 0.5, 0.0, 1.0)
            n['basal_like_program'] = max(n['basal_like_program'] - 0.01, 0.0)

        # ─── NFκB ───
        nfkb_in = n['IL6_signal'] * 0.3 + n['TGFB_signal'] * 0.15 + n['ROS_level'] * 0.2
        n['NFKB_active'] = np.clip(nfkb_in, 0.0, 1.0)
        if n['NFKB_active'] > 0.3:
            n['CXCL12_secretion'] = min(n['CXCL12_secretion'] + n['NFKB_active'] * 0.02, 1.0)
            n['CD47_expression']  = min(n['CD47_expression']  + n['NFKB_active'] * 0.01, 1.0)

        # ─── JAK/STAT3 ───
        stat3_in = n['IL6_signal'] * 0.5 + n['EGFR_active'] * 0.30 + n['HER2_active'] * 0.10
        n['STAT3_active'] = np.clip(stat3_in, 0.0, 1.0)
        if n['STAT3_active'] > 0.3:
            n['BCL2_active']     = min(n['BCL2_active']     + 0.04, 1.0)
            n['PDL1_expression'] = min(n['PDL1_expression'] + 0.02, 1.0)
            n['IDO1_active']     = min(n['IDO1_active']     + n['STAT3_active'] * 0.02, 1.0)

        # ─── MYC ───
        myc_in = (n['KRAS_active'] * 0.30 + n['ERK_active'] * 0.30
                  + n['STAT3_active'] * 0.15 + n['BCAT_nuclear'] * 0.15
                  + n['mTOR_active'] * 0.10)
        if self.mutations.get('MYC') != 'AMP':
            n['MYC_active'] = np.clip(myc_in, 0.0, 1.0)

        # ─── Metabolismo expandido ───
        n['LDHA_active'] = np.clip(
            0.15 + n['HIF1A_active'] * 0.40 + n['MYC_active'] * 0.20 + n['KRAS_active'] * 0.10,
            0.0, 1.0)
        n['warburg_effect'] = np.clip(n['warburg_effect'] * 0.7 + n['LDHA_active'] * 0.3, 0.0, 1.0)
        gls_boost = 0.30 if variant == 'G12R' else 0.10
        n['GLS_active'] = np.clip(
            0.10 + n['MYC_active'] * 0.30 + n['mTOR_active'] * 0.20 + kras_total * gls_boost,
            0.0, 1.0)
        n['FASN_active'] = np.clip(
            0.10 + n['mTOR_active'] * 0.25 + n['AKT_active'] * 0.20 + n['MYC_active'] * 0.15,
            0.0, 1.0)
        n['lipogenesis'] = np.clip(
            n['mTOR_active'] * 0.3 + n['AKT_active'] * 0.2 + n['FASN_active'] * 0.2, 0.0, 1.0)
        n['OXPHOS_active'] = np.clip(
            0.15 + (1.0 - n['warburg_effect']) * 0.30 + n['AMPK_active'] * 0.10
            - n['HIF1A_active'] * 0.10, 0.0, 1.0)
        if n['KRAS_active'] > 0.6:
            n['macropinocytosis'] = max(n['macropinocytosis'], kp['macro_boost'])
            n['autophagy']        = max(n['autophagy'], 0.3)
        if oxy < 0.3:
            n['autophagy'] = min(n['autophagy'] + 0.10, 1.0)
        if n['mTOR_active'] > 0.5:
            n['autophagy'] = max(n['autophagy'] - 0.05, 0.10)

        # ─── NRF2/KEAP1 (nuevo v9) ───
        nrf2_in = (0.10 + kras_total * kp['nrf2_boost']
                   + n['ROS_level'] * 0.15 + n['autophagy'] * 0.05)
        if self.mutations.get('ARID1A', 'WT') == 'WT':
            n['NRF2_active'] = np.clip(nrf2_in, 0.0, 1.0)

        # ─── FERROPTOSIS (nuevo v9) ───
        tp53_repression = n['TP53_active'] * 0.25 if self.mutations.get('TP53', 'WT') == 'WT' else 0.0
        slc7a11_in = (0.10 + n['NRF2_active'] * 0.45 + n['KRAS_active'] * 0.10 - tp53_repression)
        n['SLC7A11_expression'] = np.clip(slc7a11_in, 0.0, 1.0)
        n['GSH_level'] = np.clip(
            0.10 + n['SLC7A11_expression'] * 0.55 + n['GLS_active'] * 0.15, 0.0, 1.0)
        n['GPX4_level'] = np.clip(0.20 + n['GSH_level'] * 0.55, 0.0, 1.0)
        n['FSP1_level'] = np.clip(0.10 + n['NRF2_active'] * 0.20 + n['MYC_active'] * 0.10, 0.0, 1.0)
        lipid_ros_gen = 0.05 + n['ROS_level'] * 0.35 + (1.0 - oxy) * 0.15
        lipid_ros_sup = n['GPX4_level'] * 0.40 + n['FSP1_level'] * 0.20
        n['lipid_ROS'] = np.clip(lipid_ros_gen - lipid_ros_sup, 0.0, 1.0)
        gpx4_fsp1_avg = (n['GPX4_level'] + n['FSP1_level']) / 2.0
        ferrop_risk = max(0.0, n['lipid_ROS'] - gpx4_fsp1_avg * 0.8)
        n['ferroptosis_primed'] = np.clip(ferrop_risk * 2.5, 0.0, 1.0)

        # ─── ROS general ───
        n['ROS_level'] = np.clip(
            0.1 + (1.0 - oxy) * 0.25 + n['warburg_effect'] * 0.15, 0.0, 1.0)
        if n['ROS_level'] > 0.5:
            n['DDR_active'] = min(n['DDR_active'] + 0.05, 1.0)

        # ─── DDR detallado (nuevo v9) ───
        dna_damage = n['ROS_level'] * 0.3 + n['DDR_active'] * 0.5
        n['ATM_active']   = np.clip(dna_damage * 0.8, 0.0, 1.0)
        replic_stress = n['warburg_effect'] * 0.2 + n['ROS_level'] * 0.2 + n['DDR_active'] * 0.4
        n['ATR_active']   = np.clip(replic_stress * 0.75, 0.0, 1.0)
        n['CHEK1_active'] = np.clip(n['ATR_active'] * 0.70, 0.0, 1.0)
        n['CHEK2_active'] = np.clip(n['ATM_active'] * 0.70, 0.0, 1.0)
        if n['CHEK2_active'] > 0.5 and self.mutations.get('TP53') not in ('HOM_LOSS', 'HET_LOSS'):
            n['TP53_active'] = min(n['TP53_active'] + n['CHEK2_active'] * 0.08, 1.0)
        n['RAD51_active'] = np.clip(
            n['BRCA_functional'] * 0.50 + n['STAT3_active'] * 0.15 + 0.05, 0.0, 1.0)

        # ─── Inmunoevasión ampliada ───
        if n['ERK_active'] > 0.5:
            n['PDL1_expression'] = min(n['PDL1_expression'] + n['ERK_active'] * 0.03, 1.0)
        ido1_in = (0.03 + n['NFKB_active'] * 0.20 + n['ERK_active'] * 0.10 + ifng * 0.30)
        n['IDO1_active'] = np.clip(ido1_in, 0.0, 1.0)
        cd73_in = (0.08 + n['HIF1A_active'] * 0.25 + n['TGFB_signal'] * 0.20 + n['NFKB_active'] * 0.15)
        n['CD73_expression'] = np.clip(cd73_in, 0.0, 1.0)
        cd47_in = (0.15 + n['KRAS_active'] * 0.20 + n['MYC_active'] * 0.15 + n['NFKB_active'] * 0.10)
        n['CD47_expression'] = np.clip(cd47_in, 0.0, 1.0)

        # ─── SASP → NFκB amplification ───
        if n.get('SASP_active', 0) > 0.5:
            n['NFKB_active'] = min(n['NFKB_active'] + 0.15, 1.0)
            n['BCL2_active'] = min(n['BCL2_active'] + 0.10, 1.0)

        # ════════════════════════════════════════════
        # SEÑALES DE SALIDA
        # ════════════════════════════════════════════
        def hill(x, k=0.5, n_h=3):
            xc = max(x, 0.0)
            return (xc**n_h) / (k**n_h + xc**n_h + 1e-12)

        # G1→S
        g1_s = 0.05
        kras_erk = n['KRAS_active'] * n['ERK_active']
        g1_s += 0.18 * hill(kras_erk,           k=0.15, n_h=2)
        g1_s += 0.08 * hill(n['mTOR_active'],    k=0.25, n_h=2)
        g1_s += 0.04 * hill(n['AKT_active'],     k=0.30, n_h=2)
        g1_s += 0.06 * hill(n['YAP_nuclear'],    k=0.35, n_h=2)
        g1_s += 0.08 * hill(n['MYC_active'],     k=0.35, n_h=2)
        g1_s += 0.04 * hill(n['NOTCH_active'],   k=0.40, n_h=2)
        g1_s += 0.03 * hill(n['BCAT_nuclear'],   k=0.30, n_h=2)
        g1_s += 0.04 * hill(n['STAT3_active'],   k=0.25, n_h=2)
        g1_s += 0.04 * hill(n['EGFR_active'],    k=0.50, n_h=2)
        g1_s += 0.03 * hill(n['HER2_active'],    k=0.40, n_h=2)
        g1_s *= 1.0 - 0.50 * hill(n['CDKN2A_active'], k=0.40, n_h=3)
        g1_s *= 1.0 - 0.40 * hill(n['P21_active'],    k=0.40, n_h=3)
        g1_s *= 1.0 - 0.80 * hill(n['RB_active'],     k=0.40, n_h=3)
        g1_s *= 1.0 - 0.30 * hill(n['CHEK1_active'],  k=0.50, n_h=2)
        if oxy < 0.15 and n['warburg_effect'] < 0.6:
            g1_s *= 0.5
        if glucose < 0.2 and n['autophagy'] < 0.6:
            g1_s *= 0.6
        n['G1_S_transition'] = np.clip(g1_s * self.noise_prolif, 0.01, 0.90)

        # S phase
        s_prog = 0.50 + 0.2 * n['warburg_effect']
        s_prog *= 1.0 - 0.60 * hill(n['DDR_active'],   k=0.40, n_h=3)
        s_prog *= 1.0 - 0.50 * hill(n['ROS_level'],    k=0.40, n_h=3)
        s_prog *= 1.0 - 0.40 * hill(n['CHEK1_active'], k=0.35, n_h=2)
        n['S_progression'] = np.clip(s_prog * self.noise_prolif, 0.01, 0.95)

        # G2/M
        g2_m = 0.60
        g2_m *= 1.0 - 0.70 * hill(n['P21_active'],   k=0.30, n_h=3)
        g2_m *= 1.0 - 0.60 * hill(n['DDR_active'],   k=0.40, n_h=3)
        g2_m *= 1.0 - 0.40 * hill(n['CHEK2_active'], k=0.40, n_h=2)
        n['G2_M_transition'] = np.clip(g2_m * self.noise_prolif, 0.01, 0.95)

        prol_avg = (n['G1_S_transition'] + n['S_progression'] + n['G2_M_transition']) / 3.0
        n['proliferation_signal'] = np.clip(prol_avg, 0.02, 0.80)

        # Supervivencia
        surv = 0.15
        surv += 0.12 * hill(n['AKT_active'],       k=0.25, n_h=2)
        surv += 0.10 * hill(n['BCL2_active'],       k=0.40, n_h=2)
        surv += 0.06 * hill(n['autophagy'],         k=0.25, n_h=2)
        surv += 0.05 * hill(n['macropinocytosis'],  k=0.20, n_h=2)
        surv += 0.05 * hill(n['NFKB_active'],       k=0.20, n_h=2)
        surv += 0.06 * hill(n['STAT3_active'],      k=0.20, n_h=2)
        surv += 0.06 * hill(n['OXPHOS_active'],     k=0.20, n_h=2)
        surv += 0.04 * hill(n['DDR_active'],        k=0.30, n_h=2)
        surv += 0.05 * hill(n['TEAD_active'],       k=0.35, n_h=2)
        surv += 0.10 * hill(n['EGFR_active'],       k=0.40, n_h=2)
        surv += 0.08 * hill(n['HER2_active'],       k=0.35, n_h=2)
        surv += 0.05 * hill(n['GPX4_level'],        k=0.40, n_h=2)
        n['survival_signal'] = np.clip(surv * self.noise_surv, 0.1, 0.85)

        # Apoptosis
        apop = 0.002
        if n['TP53_active'] > 0.5:
            stress = max(n['DDR_active'], n['ROS_level'], 1.0 - oxy)
            if stress > 0.3:
                apop += 0.02 * (stress - 0.3)
        if oxy < 0.1:
            apop += 0.02
        bax_bcl2 = n['BAX_active'] - n['BCL2_active']
        if bax_bcl2 > 0:
            apop += bax_bcl2 * 0.05
        if n['DDR_active'] > 0.5 and n['BRCA_functional'] < 0.3:
            apop += 0.04
        if n['BRCA_functional'] < 0.3 and n['ROS_level'] > 0.35:
            apop += (n['ROS_level'] - 0.35) * 0.10
        if n['ATM_active'] > 0.6 and n['BRCA_functional'] < 0.2:
            apop += 0.03
        # KRASi withdrawal shock
        if variant != 'WT' and self.mutations.get('KRAS') == 'HET_MUT':
            kras_mutant = n['KRAS_active']
            kras_baseline = self.kras_profile['kras_level']
            kras_suppression = max(0, 1.0 - (kras_mutant / max(kras_baseline, 0.01)))
            if kras_suppression > 0.3:
                shock_decay = max(0.0, 1.0 - temporal_factor * 0.8)
                apop += kras_suppression * 0.20 * shock_decay
        if n['ROS_level'] > 0.6:
            apop += 0.01
        apop *= (1.0 - n['survival_signal'] * 0.3)
        n['apoptosis_signal'] = np.clip(apop * self.noise_apop, 0.001, 0.20)

    # ─── Getters ───
    def get_proliferation_rate(self):
        return self.nodes['proliferation_signal']

    def get_apoptosis_probability(self):
        return self.nodes['apoptosis_signal']

    def get_survival(self):
        return self.nodes['survival_signal']

    def get_ferroptosis_risk(self):
        return self.nodes.get('ferroptosis_primed', 0.0)

    def is_basal_like(self):
        return self.nodes['basal_like_program'] > 0.5

    def get_pdl1_level(self):
        return self.nodes['PDL1_expression']

    def get_state_summary(self):
        n = self.nodes
        return {
            'KRAS': n['KRAS_active'], 'RAF': n['RAF_active'],
            'MEK': n['MEK_active'], 'ERK': n['ERK_active'],
            'PI3K': n['PI3K_active'], 'AKT': n['AKT_active'],
            'mTOR': n['mTOR_active'], 'YAP': n['YAP_nuclear'],
            'TEAD': n['TEAD_active'], 'HIF-1α': n['HIF1A_active'],
            'HIF-2α': n['HIF2A_active'], 'STAT3': n['STAT3_active'],
            'NFκB': n['NFKB_active'], 'WNT': n['WNT_active'],
            'p53': n['TP53_active'], 'p21': n['P21_active'],
            'PD-L1': n['PDL1_expression'], 'MHC-I': n['MHC1_expression'],
            'E-cad': n['ECAD_expression'], 'Vim': n['VIM_expression'],
            'Autofagia': n['autophagy'], 'Warburg': n['warburg_effect'],
            'Macropino': n['macropinocytosis'],
            'ROS': n['ROS_level'], 'BCL2': n['BCL2_active'],
            'MYC': n['MYC_active'], 'OXPHOS': n['OXPHOS_active'],
            'EGFR': n['EGFR_active'], 'HER2': n['HER2_active'],
            'NRF2': n['NRF2_active'],
            'GPX4': n['GPX4_level'], 'SLC7A11': n['SLC7A11_expression'],
            'GSH': n['GSH_level'], 'lipid_ROS': n['lipid_ROS'],
            'Ferroptosis_risk': n['ferroptosis_primed'],
            'AMPK': n['AMPK_active'], 'LKB1': n['LKB1_active'],
            'LDHA': n['LDHA_active'], 'GLS': n['GLS_active'],
            'IDO1': n['IDO1_active'], 'CD73': n['CD73_expression'],
            'CD47': n['CD47_expression'],
            'ATM': n['ATM_active'], 'ATR': n['ATR_active'],
            'CHK1': n['CHEK1_active'], 'RAD51': n['RAD51_active'],
            'Prolif': n['proliferation_signal'],
            'Superviv': n['survival_signal'],
            'Apoptosis': n['apoptosis_signal'],
            'KRAS_variant': self.kras_variant,
            'Subtipo': 'basal-like' if self.is_basal_like() else 'classical',
        }


# ════════════════════════════════════════════════════════
# GENERADOR DE PERFIL MUTACIONAL (frecuencias TCGA/ICGC)
# ════════════════════════════════════════════════════════
def generate_mutation_profile(rng=None) -> dict:
    """
    Perfil mutacional PDAC v9. Frecuencias ajustadas a TCGA (n=185) + ICGC PACA-AU (n=396).
    Incluye LKB1/STK11, ARID1A.
    """
    if rng is None:
        rng = np.random.default_rng()

    has_kras = rng.random() < 0.92
    variant = 'WT'
    if has_kras:
        r = rng.random()
        if r < 0.41:    variant = 'G12D'
        elif r < 0.71:  variant = 'G12V'
        elif r < 0.87:  variant = 'G12R'
        elif r < 0.89:  variant = 'G12C'
        elif r < 0.91:  variant = 'Q61H'
        else:           variant = 'G12D'
        kras_status = 'HET_MUT'
    else:
        kras_status = 'WT'

    def get_loss_status(prob_loss, prob_hom_given_loss):
        if rng.random() < prob_loss:
            return 'HOM_LOSS' if rng.random() < prob_hom_given_loss else 'HET_LOSS'
        return 'WT'

    return {
        'KRAS':       kras_status,
        'KRAS_variant': variant,
        'TP53':       get_loss_status(0.70, 0.60),
        'CDKN2A':     get_loss_status(0.52, 0.70),
        'SMAD4':      get_loss_status(0.55, 0.50),
        'YAP_amp':    'AMP' if rng.random() < 0.30 else 'WT',
        'BRCA':       get_loss_status(0.07, 0.70),
        'PTEN':       get_loss_status(0.05, 0.50),
        'MYC':        'AMP' if rng.random() < 0.15 else 'WT',
        'LKB1':       get_loss_status(0.05, 0.60),
        'ARID1A':     get_loss_status(0.07, 0.40),
        'MSI_status': 'MSI-H' if rng.random() < 0.02 else 'MSS',
    }
