"""
Agente: Célula cancerosa PDAC v7.
- Heterogeneidad clonal (mutaciones per-cell + evolución somática)
- Resistencia adquirida mutacional
- Farmacocinética (penetración ECM)
Mesa 3.x compatible.
"""

import mesa
import numpy as np
from signaling.pdac_network import SignalingNetwork, generate_mutation_profile


class CancerCell(mesa.Agent):
    """Célula cancerosa con señalización molecular y heterogeneidad clonal."""

    def __init__(self, model, _pos=None, mutations=None, generation=0):
        super().__init__(model)
        self.cell_type = 'cancer'

        # Perfil mutacional con heterogeneidad clonal
        if mutations is None:
            mutations = generate_mutation_profile(rng=np.random.default_rng())
        self.mutations = dict(mutations)  # Copia propia
        self.generation = generation  # Generación clonal
        self.alive = True
        self.subtype = 'classical'
        self.divisions = 0
        self.age = 0
        self.drug_exposure = {}
        self.drug_exposure_total = 0.0
        self.resistant = False
        self.resistance_mutations = []

        self.cell_cycle_phase = 'G1'
        # Desincronizar ciclo celular inicial para evitar divisiones explosivas simultáneas
        self.cycle_timer = float(model.random.uniform(0.0, 10.0)) if hasattr(model, 'random') else 0.0

        # Estado de senescencia (TIS o OIS)
        self.senescent = False
        self.senescence_counter = 0      # Steps acumulados en estado senescente
        self.p21_accumulation = 0.0      # Proxy de p21 por daño sublethal (→ TIS)

        # Variantes subclonales (cada célula puede tener alteraciones extra)
        if generation > 0 and hasattr(model, '_rng'):
            self._apply_somatic_evolution(model._rng)

        self.signaling = SignalingNetwork(mutations=self.mutations)

    def _apply_somatic_evolution(self, rng):
        """
        Evolución somática: cada división puede generar alteraciones nuevas.
        Tasa ~1 mutación por 10⁸ bases / división → simplificado a probabilidades
        por vía oncogénica.
        """
        # Probabilidad baja de ganar mutaciones driver nuevas por generación
        gen_adj = min(self.generation, 20)
        if rng.random() < 0.005 * gen_adj:
            # Puede ganar amplificación MYC
            if self.mutations.get('MYC') != 'AMP' and rng.random() < 0.02:
                self.mutations['MYC'] = 'AMP'
                self.resistance_mutations.append('MYC_amplification')
        if rng.random() < 0.003 * gen_adj:
            if self.mutations.get('YAP_amp') != 'AMP' and rng.random() < 0.03:
                self.mutations['YAP_amp'] = 'AMP'
                self.resistance_mutations.append('YAP_amplification')
        if rng.random() < 0.002 * gen_adj:
            if self.mutations.get('PTEN') not in ('HOM_LOSS', 'HET_LOSS') and rng.random() < 0.02:
                self.mutations['PTEN'] = 'HET_LOSS'
                self.resistance_mutations.append('PTEN_loss')
        # Mutaciones de resistencia a KRASi (Cancer Discovery 2024)
        if rng.random() < 0.001 * gen_adj:
            if rng.random() < 0.03:
                self.resistance_mutations.append('NRAS_Q61K')
        if rng.random() < 0.002 * gen_adj:
            if rng.random() < 0.02:
                self.resistance_mutations.append('EGFR_amplification')

    def step(self):
        if not self.alive:
            return
        self.age += 1

        # Microambiente local
        microenv = {}
        if hasattr(self.model, 'microenv') and self.pos is not None:
            microenv = self.model.microenv.get_local_conditions(
                self.pos[0], self.pos[1])

        # Fármacos — con penetración ECM
        phases_fx = getattr(self.model, 'drug_effects_by_phase', None)
        if phases_fx:
            drug_fx = phases_fx.get(self.cell_cycle_phase, {}).copy()
        else:
            drug_fx = getattr(self.model, 'drug_effects', {}).copy()

        if drug_fx and hasattr(self.model, 'microenv') and self.pos is not None:
            pen = self.model.microenv.drug_penetration_factor(
                self.pos[0], self.pos[1])
            if pen < 1.0:
                drug_fx = {k: v * pen for k, v in drug_fx.items()}

        # Farmacocinética real (µM) ya está incorporada en current_doses_uM.
        # Fallback solo por si acaso:
        if drug_fx and hasattr(self.model, 'pk_factor') and not hasattr(self.model, 'current_doses_uM'):
            pk = self.model.pk_factor
            drug_fx = {k: v * pk for k, v in drug_fx.items()}

        # Resistencia adquirida reduce efecto de fármacos
        if self.resistant and drug_fx:
            drug_fx = {k: v * 0.4 for k, v in drug_fx.items()}  # 60% resistencia

        # Señalización
        active_drugs = getattr(self.model, 'active_drugs', [])
        # Pasamos cell_cycle_phase explicitamente a la red
        self.signaling.update(microenv=microenv, drug_effects=drug_fx, active_drugs=active_drugs, cell_cycle_phase=self.cell_cycle_phase)
        self.subtype = 'basal-like' if self.signaling.is_basal_like() else 'classical'

        rng = self.random

        # ── Apoptosis ──────────────────────────────────────────────────────────
        apop_prob = self.signaling.get_apoptosis_probability()
        # Células resistentes tienen menos apoptosis
        if self.resistant:
            apop_prob *= 0.5
        # Células senescentes: BCL-2 las protege de apoptosis basal,
        # pero navitoclax las hace hipersensibles (senolisis)
        if self.senescent:
            apop_prob *= 0.3   # Resistentes a apoptosis basal
            navitoclax_fx = drug_fx.get('SASP_active', 0.0)   # Navitoclax target SASP_active
            bcl2_inh = drug_fx.get('BCL2_active', 0.0)
            if bcl2_inh > 0.3:
                # Navitoclax → senolisis. Células senescentes mueren preferentemente.
                senolysis_p = bcl2_inh * 0.70
                if rng.random() < senolysis_p:
                    self._die()
                    return
        if rng.random() < apop_prob:
            self._die()
            return

        # ── Ferroptosis ────────────────────────────────────────────────────────
        # Muerte por lipid ROS cuando GPX4 y FSP1 son insuficientes
        ferrop_risk = self.signaling.get_ferroptosis_risk()
        if ferrop_risk > 0.55:
            ferrop_prob = (ferrop_risk - 0.55) * 0.18   # Prob creciente con el riesgo
            if rng.random() < ferrop_prob:
                self._die()
                return

        # ── Senescencia (TIS / OIS) ────────────────────────────────────────────
        if not self.senescent:
            # TIS: daño sublethal acumulado por quimio → p21 → G1 arrest con SASP.
            # La condición captura tanto DDR de red de señalización como señal directa de fármacos.
            ddr_signal = self.signaling.nodes.get('DDR_active', 0.0)
            atr_signal = self.signaling.nodes.get('ATR_active', 0.0)
            survival   = self.signaling.nodes.get('survival_signal', 0.5)
            # Drug-induced DNA damage: apoptosis_signal<0 o DDR_active<0 en drug_fx significa activación
            drug_dna_dmg = abs(min(drug_fx.get('apoptosis_signal', 0.0), 0.0))
            drug_ddr_act = abs(min(drug_fx.get('DDR_active', 0.0), 0.0))
            drug_atr_act = abs(min(drug_fx.get('ATR_active', 0.0), 0.0))
            # DDR efectivo = mayor entre señal de red y señal directa de fármaco
            effective_ddr = max(ddr_signal, drug_ddr_act, drug_atr_act, drug_dna_dmg * 0.7)
            if effective_ddr > 0.15 and 0.10 < survival < 0.65 and drug_fx:
                self.p21_accumulation += effective_ddr * 0.08
            # Umbral de senescencia
            if self.p21_accumulation > 0.45 and rng.random() < 0.10:
                self.senescent = True
                self.cell_cycle_phase = 'G1'   # Arresto permanente en G1
                self.cycle_timer = 0.0
        else:
            # Célula senescente: secreta SASP (IL-6 → NFκB en vecinas)
            self.senescence_counter += 1
            self._secrete_sasp()
            # Senescencia muy prolongada → eventual muerte (10-30% por step tardío)
            if self.senescence_counter > 80 and rng.random() < 0.003:
                self._die()
                return

        # Células senescentes no proliferan (bloqueo permanente G1/p21)
        if self.senescent:
            self._move()
            return

        # ── Proliferación — Ciclo Celular Estructurado (G1 -> S -> G2 -> M) ──
        growth_inh = getattr(self.model, 'growth_inhibition', 0.0)
        density_brake = (1.0 - growth_inh)  # Arresto por contacto afecta fases

        if self.cell_cycle_phase == 'G1':
            rate = self.signaling.nodes.get('G1_S_transition', 0.0) * density_brake
            self.cycle_timer += rate
            if self.cycle_timer >= 12.0:  # ~12 horas en G1
                self.cell_cycle_phase = 'S'
                self.cycle_timer = 0.0
                
        elif self.cell_cycle_phase == 'S':
            rate = self.signaling.nodes.get('S_progression', 0.0) * density_brake
            self.cycle_timer += rate
            if self.cycle_timer >= 8.0:   # ~8 horas síntesis ADN
                self.cell_cycle_phase = 'G2'
                self.cycle_timer = 0.0
                
        elif self.cell_cycle_phase == 'G2':
            rate = self.signaling.nodes.get('G2_M_transition', 0.0) * density_brake
            self.cycle_timer += rate
            if self.cycle_timer >= 4.0:   # ~4 horas preparación división
                self.cell_cycle_phase = 'M'
                self.cycle_timer = 0.0
                
        elif self.cell_cycle_phase == 'M':
            rate = 0.8  # Mitosis rápida (~1-2h), constante
            self.cycle_timer += rate
            if self.cycle_timer >= 2.0:
                success = self._proliferate()
                if success:
                    # Tras división exitosa vuelve a G1 (hija también nace en G1)
                    self.cell_cycle_phase = 'G1'
                    self.cycle_timer = 0.0
                # Si no hay espacio, queda bloqueada en catástrofe/arresto mitótico

        # Movimiento
        self._move()

        # Tracking de exposición a fármacos
        for dn in getattr(self.model, 'active_drugs', []):
            dose = self.model.drug_doses.get(dn, 0)
            if dose > 0:
                self.drug_exposure[dn] = self.drug_exposure.get(dn, 0) + 1
                self.drug_exposure_total += dose

        # Adquisición de resistencia (por exposición prolongada)
        self._check_acquired_resistance()

    def _proliferate(self):
        if self.pos is None:
            return False
        empty = self._find_empty()
        if empty is None:
            return False

        # Hija hereda mutaciones + posible evolución somática
        daughter_muts = self.mutations.copy()
        daughter = CancerCell(self.model, _pos=empty,
                              mutations=daughter_muts,
                              generation=self.generation + 1)
        # Hereda parte de la resistencia
        daughter.drug_exposure = {k: int(v * 0.8) for k, v in self.drug_exposure.items()}
        daughter.resistant = self.resistant and self.random.random() < 0.95
        daughter.resistance_mutations = list(self.resistance_mutations)
        # Hereda parcialmente el ruido epigenético (50% madre + 50% nuevo)
        daughter.signaling.noise_prolif = 0.5 * self.signaling.noise_prolif + 0.5 * daughter.signaling.noise_prolif
        daughter.signaling.noise_surv = 0.5 * self.signaling.noise_surv + 0.5 * daughter.signaling.noise_surv
        daughter.signaling.noise_apop = 0.5 * self.signaling.noise_apop + 0.5 * daughter.signaling.noise_apop
        self.model.grid.place_agent(daughter, empty)
        self.divisions += 1
        return True

    def _secrete_sasp(self):
        """
        SASP (Senescence-Associated Secretory Phenotype): secreción de IL-6, MMPs.
        Activa NFκB en células vecinas → puede promover supervivencia de células no senescentes
        (paradoja de SASP: pro-tumoral en contexto de quimio).
        """
        if self.pos is None or not hasattr(self.model, 'microenv'):
            return
        # Amplitud SASP proporcional a senescence_counter (más SASP cuanto más tiempo senescente)
        sasp_strength = min(0.15 + self.senescence_counter * 0.003, 0.55)
        x, y = self.pos
        # Secretar IL-6 al microambiente (modelado como IL6 field si existe, o como boost NFκB)
        if hasattr(self.model.microenv, 'il6') and hasattr(self.model.microenv.il6, '__setitem__'):
            try:
                self.model.microenv.il6[x, y] = min(
                    self.model.microenv.il6[x, y] + sasp_strength * 0.10, 1.0
                )
            except Exception:
                pass
        # Actualizar SASP_active en la propia red de señalización
        self.signaling.nodes['SASP_active'] = min(sasp_strength, 1.0)

    def _check_acquired_resistance(self):
        """
        Resistencia adquirida basada en duración de exposición + presión selectiva.
        Mecanismos modelados:
        - Upregulation MDR1/ABCB1 (efflux de quimio)
        - Amplificación diana (ej: KRAS Y96D post-sotorasib)
        - Activación de vías bypass
        """
        if self.resistant:
            return

        # Tras exposición (>10h), probabilidad creciente de resistencia adquirida
        for drug, hours in self.drug_exposure.items():
            if hours > 10:
                p_resist = min(0.0015 * (hours - 10), 0.06)
                if self.random.random() < p_resist:
                    self.resistant = True
                    # Tipo de resistencia según fármaco
                    if drug in ('gemcitabine', '5fu', 'oxaliplatin', 'irinotecan'):
                        self.resistance_mutations.append(f'MDR1_upregulation_{drug}')
                    elif drug in ('daraxonrasib', 'mrtx1133', 'sotorasib', 'rmc7977'):
                        self.resistance_mutations.append(f'KRAS_secondary_mutation_{drug}')
                        # Bypass vía YAP
                        self.mutations['YAP_amp'] = 'AMP'
                        self.signaling = SignalingNetwork(mutations=self.mutations)
                    elif drug == 'trametinib':
                        self.resistance_mutations.append('MEK_amplification')
                    elif drug in ('anti_pd1', 'anti_ctla4'):
                        self.resistance_mutations.append('B2M_loss_immune_escape')
                    elif drug in ('rmc4630', 'bi3406'):
                        self.resistance_mutations.append(f'RAS_amplification_{drug}')
                    elif drug == 'palbociclib':
                        self.resistance_mutations.append('CDK2_bypass')
                    elif drug == 'ceralasertib':
                        self.resistance_mutations.append('CHK1_upregulation')
                        self.signaling.nodes['CHEK1_active'] = min(
                            self.signaling.nodes.get('CHEK1_active', 0.1) + 0.30, 1.0)
                    elif drug == 'rsl3':
                        # GPX4i resistencia: NRF2 upregulation → más GPX4/FSP1
                        self.resistance_mutations.append('NRF2_upregulation_GPX4i')
                        self.mutations['NRF2'] = 'GAIN'
                        self.signaling = SignalingNetwork(mutations=self.mutations)
                    elif drug in ('erastin', 'sulfasalazine'):
                        # xCT/SLC7A11i resistencia: transcripción NRF2/ARE → SLC7A11
                        self.resistance_mutations.append('SLC7A11_upregulation_xCTi')
                        self.signaling.nodes['NRF2_active'] = min(
                            self.signaling.nodes.get('NRF2_active', 0.2) + 0.35, 1.0)
                    elif drug == 'navitoclax':
                        # BCL-2i resistencia: upregulation MCL-1 (modelado via survival + NFκB)
                        self.resistance_mutations.append('MCL1_upregulation')
                        self.signaling.nodes['survival_signal'] = min(
                            self.signaling.nodes.get('survival_signal', 0.3) + 0.25, 1.0)
                    elif drug == 'metformin':
                        self.resistance_mutations.append('OXPHOS_bypass_metformin')
                    break

    def _move(self):
        if self.pos is None:
            return
        if hasattr(self.model, 'microenv'):
            if self.model.microenv.ecm_blocks_movement(self.pos[0], self.pos[1]):
                return
        # EMT → más móviles
        move_prob = 0.15
        if self.signaling.nodes.get('invasion_signal', 0) > 0.3:
            move_prob = 0.35  # Células mesenquimales más móviles
        if self.random.random() < move_prob:
            empty = self._find_empty()
            if empty is not None:
                self.model.grid.move_agent(self, empty)

    def _die(self):
        self.alive = False
        if self.pos is not None:
            self.model.grid.remove_agent(self)

    def _find_empty(self):
        if self.pos is None:
            return None
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        empty = [p for p in nbrs if self.model.grid.is_cell_empty(p)]
        return self.random.choice(empty) if empty else None
