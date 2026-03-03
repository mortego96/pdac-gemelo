"""
Agentes inmunes PDAC v7.
- CD8+ T-cell con exhaustion progresivo (Tim-3, LAG-3, PD-1)
- Macrófago M2 (protumoral)
- Treg (CD4+FOXP3+, supresor)
- MDSC (mieloid-derived suppressor cell)
- NK cell (natural killer)
Mesa 3.x compatible.
"""

import mesa


class CD8_Tcell(mesa.Agent):
    """
    Linfocito T CD8+ citotóxico.
    Modelo de exhaustion progresivo:
    - Fase 1 (competente): mata células cancerosas con PD-L1 bajo
    - Fase 2 (pre-exhausto): PD-1+, eficacia reducida
    - Fase 3 (exhausto): PD-1+ Tim-3+ LAG-3+, casi inactivo
    Anti-PD-1 puede revertir parcialmente la exhaustion.
    (Wherry EJ, Nat Immunol 2011; Philip M, Science 2017)
    """

    def __init__(self, model, pos=None):
        super().__init__(model)
        self.cell_type = 'tcell'
        self.alive = True
        self.kills = 0
        self.exhaustion = 0.0      # 0=competente, 1=exhausto
        self.pd1_level = 0.1       # Expresión PD-1
        self.tim3_level = 0.0      # Expresión Tim-3
        self.lag3_level = 0.0      # Expresión LAG-3
        self.ifng_production = 1.0 # Capacidad de producir IFN-γ
        self.activated = False
        self.age = 0

    def step(self):
        if not self.alive:
            return
        self.age += 1

        me = {}
        if hasattr(self.model, 'microenv') and self.pos is not None:
            me = self.model.microenv.get_local_conditions(self.pos[0], self.pos[1])

        # Exhaustion progresivo por exposición al TME
        self._update_exhaustion(me)

        # Anti-PD-1 revierte parcialmente la exhaustion, pero mucho menos si es MSS
        anti_pd1 = self.model.drug_doses.get('anti_pd1', 0)
        msi = getattr(self.model, 'msi_status', 'MSS')
        if msi == 'MSS':
            anti_pd1 *= 0.1  # 90% reducción de eficacia en pMMR/MSS (gran limitación PDAC)

        if anti_pd1 > 0 and self.exhaustion < 0.8:
            self.pd1_level *= (1.0 - anti_pd1 * 0.6)
            self.exhaustion = max(self.exhaustion - anti_pd1 * 0.1, 0.0)
            self.ifng_production = min(self.ifng_production + 0.05, 1.0)

        # Anti-CTLA-4
        anti_ctla4 = self.model.drug_doses.get('anti_ctla4', 0)
        if msi == 'MSS':
            anti_ctla4 *= 0.1

        if anti_ctla4 > 0:
            self.activated = True
            self.ifng_production = min(self.ifng_production + 0.03 * anti_ctla4, 1.0)

        # Intentar matar
        efficacy = self._killing_efficacy()
        if efficacy > 0.05:
            self._attempt_kill(efficacy)

        # Secretar IFN-γ
        if self.activated and self.pos is not None and hasattr(self.model, 'microenv'):
            me_obj = self.model.microenv
            x, y = self.pos
            if 0 <= x < me_obj.width and 0 <= y < me_obj.height:
                me_obj.ifng[x, y] = min(
                    me_obj.ifng[x, y] + 0.02 * self.ifng_production, 1.0)

        # Movimiento (quimiotaxis hacia tumor)
        self._move(me)

        # Muerte por TME hostil
        if me.get('oxygen', 0.8) < 0.1 or self.exhaustion > 0.95:
            if self.random.random() < 0.05:
                self._die()

    def _update_exhaustion(self, me):
        """Exhaustion progresivo por señales del TME."""
        # TGF-β, IL-6 y hipoxia → agotan a los T-cells
        tgfb = me.get('tgfb', 0)
        il6 = me.get('il6', 0)
        lactate = me.get('lactate', 0)

        exhaust_rate = 0.002  # Basal
        exhaust_rate += tgfb * 0.005
        exhaust_rate += il6 * 0.003
        exhaust_rate += lactate * 0.004  # Acidosis → supresión T-cell
        self.exhaustion = min(self.exhaustion + exhaust_rate, 1.0)

        # Actualizar marcadores de exhaustion
        if self.exhaustion > 0.3:
            self.pd1_level = min(0.2 + self.exhaustion * 0.6, 1.0)
            self.activated = True
        if self.exhaustion > 0.5:
            self.tim3_level = min((self.exhaustion - 0.5) * 1.5, 1.0)
        if self.exhaustion > 0.7:
            self.lag3_level = min((self.exhaustion - 0.7) * 2.0, 1.0)

        # IFN-γ baja con exhaustion
        self.ifng_production = max(1.0 - self.exhaustion * 0.8, 0.1)

    def _killing_efficacy(self):
        """Eficacia citotóxica basada en estado de exhaustion."""
        base = 0.3
        # Exhaustion reduce eficacia
        base *= (1.0 - self.exhaustion * 0.8)
        # Necesita estar activado
        if not self.activated:
            base *= 0.5
        return max(base, 0.02)

    def _attempt_kill(self, efficacy):
        if self.pos is None:
            return

        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        for pos in nbrs:
            agents = self.model.grid.get_cell_list_contents([pos])
            for a in agents:
                if (hasattr(a, 'cell_type') and a.cell_type == 'cancer'
                        and a.alive):
                    # PD-L1 del tumor bloquea killing
                    if hasattr(a, 'signaling'):
                        pdl1 = a.signaling.get_pdl1_level()
                        block = pdl1 * self.pd1_level
                        efficacy *= (1.0 - block * 0.7)
                        # MHC-I necesario para reconocimiento
                        mhc1 = a.signaling.nodes.get('MHC1_expression', 0.7)
                        if mhc1 < 0.3:
                            efficacy *= 0.2  # Casi invisible

                    if self.random.random() < efficacy:
                        a._die()
                        self.kills += 1
                        self.activated = True
                        # Matar agota parcialmente
                        self.exhaustion = min(self.exhaustion + 0.02, 1.0)
                        return

    def _move(self, me):
        if self.pos is None:
            return
        # ECM bloqueo
        if hasattr(self.model, 'microenv'):
            if self.model.microenv.ecm_blocks_movement(self.pos[0], self.pos[1], 0.6):
                return
        if self.random.random() < 0.2:
            nbrs = self.model.grid.get_neighborhood(
                self.pos, moore=True, include_center=False)
            empty = [p for p in nbrs if self.model.grid.is_cell_empty(p)]
            if empty:
                self.model.grid.move_agent(self, self.random.choice(empty))

    def _die(self):
        self.alive = False
        if self.pos is not None:
            self.model.grid.remove_agent(self)
        self.remove()


class MacrophageM2(mesa.Agent):
    """
    Macrófago M2 (tumor-associated macrophage, TAM).
    Protumoral: secreta IL-6, TGF-β, suprime T-cells.
    """

    def __init__(self, model, pos=None):
        super().__init__(model)
        self.cell_type = 'macrophage'
        self.alive = True
        self.polarization = 'M2'  # M1 (antitumoral) o M2 (protumoral)
        self.age = 0

    def step(self):
        if not self.alive:
            return
        self.age += 1

        me = {}
        if hasattr(self.model, 'microenv') and self.pos is not None:
            me = self.model.microenv.get_local_conditions(self.pos[0], self.pos[1])
            x, y = self.pos
            me_obj = self.model.microenv

            # IFN-γ del TME puede repolarizar a M1
            if me.get('ifng', 0) > 0.4:
                if self.random.random() < 0.03:
                    self.polarization = 'M1'
            # Lactato → M2
            if me.get('lactate', 0) > 0.4:
                self.polarization = 'M2'

            if self.polarization == 'M2':
                # Secretar IL-6 y TGF-β (protumoral)
                if 0 <= x < me_obj.width and 0 <= y < me_obj.height:
                    me_obj.il6[x, y] = min(me_obj.il6[x, y] + 0.02, 1.0)
                    me_obj.tgfb[x, y] = min(me_obj.tgfb[x, y] + 0.015, 1.0)
                    me_obj.tnfa[x, y] = min(me_obj.tnfa[x, y] + 0.01, 1.0)
            elif self.polarization == 'M1':
                # M1: fagocitosis antitumoral (limitada)
                self._attempt_phagocytosis()

        # Movimiento
        if self.random.random() < 0.1:
            self._move()

    def _attempt_phagocytosis(self):
        """M1 puede fagocitar células tumorales (baja eficacia)."""
        if self.pos is None:
            return
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        for pos in nbrs:
            agents = self.model.grid.get_cell_list_contents([pos])
            for a in agents:
                if (hasattr(a, 'cell_type') and a.cell_type == 'cancer'
                        and a.alive):
                    if self.random.random() < 0.03:  # Baja eficacia
                        a._die()
                        return

    def _move(self):
        if self.pos is None:
            return
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        empty = [p for p in nbrs if self.model.grid.is_cell_empty(p)]
        if empty:
            self.model.grid.move_agent(self, self.random.choice(empty))


class Treg(mesa.Agent):
    """
    Linfocito T regulador (CD4+FOXP3+).
    Inmunosupresor: inhibe CD8+ y NK directamente.
    Reclutado por CXCL12 y TGF-β de CAFs.
    Anti-CTLA-4 puede deplecionar Tregs parcialmente.
    """

    def __init__(self, model, pos=None):
        super().__init__(model)
        self.cell_type = 'treg'
        self.alive = True
        self.suppression_radius = 2
        self.age = 0

    def step(self):
        if not self.alive:
            return
        self.age += 1

        # Anti-CTLA-4 puede deplecionar Tregs, pero capado en MSS
        anti_ctla4 = self.model.drug_doses.get('anti_ctla4', 0)
        msi = getattr(self.model, 'msi_status', 'MSS')
        if msi == 'MSS':
            anti_ctla4 *= 0.1

        if anti_ctla4 > 0:
            if self.random.random() < anti_ctla4 * 0.05:
                self._die()
                return

        # Suprimir T-cells cercanos → aumentar su exhaustion
        self._suppress_tcells()

        # Secretar citoquinas inmunosupresoras
        if self.pos is not None and hasattr(self.model, 'microenv'):
            x, y = self.pos
            me = self.model.microenv
            if 0 <= x < me.width and 0 <= y < me.height:
                me.tgfb[x, y] = min(me.tgfb[x, y] + 0.02, 1.0)
                me.il6[x, y] = min(me.il6[x, y] + 0.01, 1.0)

        # Movimiento
        if self.random.random() < 0.08:
            self._move()

    def _suppress_tcells(self):
        """Suprime CD8+ y NK cercanos → aumenta exhaustion."""
        if self.pos is None:
            return
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False, radius=self.suppression_radius)
        for pos in nbrs:
            agents = self.model.grid.get_cell_list_contents([pos])
            for a in agents:
                if hasattr(a, 'cell_type') and a.cell_type == 'tcell' and a.alive:
                    a.exhaustion = min(a.exhaustion + 0.01, 1.0)
                elif hasattr(a, 'cell_type') and a.cell_type == 'nk' and a.alive:
                    a.suppressed = True

    def _move(self):
        if self.pos is None:
            return
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        empty = [p for p in nbrs if self.model.grid.is_cell_empty(p)]
        if empty:
            self.model.grid.move_agent(self, self.random.choice(empty))

    def _die(self):
        self.alive = False
        if self.pos is not None:
            self.model.grid.remove_agent(self)
        self.remove()


class MDSC(mesa.Agent):
    """
    Myeloid-Derived Suppressor Cell.
    Suprime inmunidad adaptativa. Frecuente en PDAC.
    Reclutada por GM-CSF y CXCL12 del tumor.
    """

    def __init__(self, model, pos=None):
        super().__init__(model)
        self.cell_type = 'mdsc'
        self.alive = True
        self.age = 0

    def step(self):
        if not self.alive:
            return
        self.age += 1

        # Suprimir T-cells cercanos
        if self.pos is not None:
            nbrs = self.model.grid.get_neighborhood(
                self.pos, moore=True, include_center=False, radius=2)
            for pos in nbrs:
                agents = self.model.grid.get_cell_list_contents([pos])
                for a in agents:
                    if hasattr(a, 'cell_type') and a.cell_type == 'tcell' and a.alive:
                        a.exhaustion = min(a.exhaustion + 0.008, 1.0)

        # Producir ROS y arginasa (suprime T-cells)
        if self.pos is not None and hasattr(self.model, 'microenv'):
            x, y = self.pos
            me = self.model.microenv
            if 0 <= x < me.width and 0 <= y < me.height:
                me.ros[x, y] = min(me.ros[x, y] + 0.01, 1.0)

        # Movimiento
        if self.random.random() < 0.1:
            self._move()

    def _move(self):
        if self.pos is None:
            return
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        empty = [p for p in nbrs if self.model.grid.is_cell_empty(p)]
        if empty:
            self.model.grid.move_agent(self, self.random.choice(empty))


class NKCell(mesa.Agent):
    """
    Natural Killer cell.
    Mata células con MHC-I bajo (complementa CD8+).
    Suprimida por TGF-β y Tregs.
    """

    def __init__(self, model, pos=None):
        super().__init__(model)
        self.cell_type = 'nk'
        self.alive = True
        self.suppressed = False
        self.kills = 0
        self.age = 0

    def step(self):
        if not self.alive:
            return
        self.age += 1
        self.suppressed = False  # Reset each step

        me = {}
        if hasattr(self.model, 'microenv') and self.pos is not None:
            me = self.model.microenv.get_local_conditions(self.pos[0], self.pos[1])

        # TGF-β inhibe NK
        if me.get('tgfb', 0) > 0.5:
            self.suppressed = True

        # Intentar matar (mejor contra MHC-I bajo — opuesto a CD8+)
        if not self.suppressed:
            self._attempt_kill()

        # Movimiento
        if self.random.random() < 0.15:
            self._move()

    def _attempt_kill(self):
        if self.pos is None:
            return
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        for pos in nbrs:
            agents = self.model.grid.get_cell_list_contents([pos])
            for a in agents:
                if (hasattr(a, 'cell_type') and a.cell_type == 'cancer'
                        and a.alive):
                    # NK mata preferentemente células con MHC-I BAJO
                    # (opuesto a CD8+ que necesita MHC-I)
                    mhc1 = 0.7
                    if hasattr(a, 'signaling'):
                        mhc1 = a.signaling.nodes.get('MHC1_expression', 0.7)
                    kill_prob = 0.1 * (1.0 - mhc1)  # Mejor si MHC-I bajo
                    if self.random.random() < kill_prob:
                        a._die()
                        self.kills += 1
                        return

    def _move(self):
        if self.pos is None:
            return
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        empty = [p for p in nbrs if self.model.grid.is_cell_empty(p)]
        if empty:
            self.model.grid.move_agent(self, self.random.choice(empty))
