"""
Agente: Fibroblastos asociados al cáncer (CAFs) v7.
3 subtipos: myCAF (αSMA+), iCAF (IL-6+), apCAF (MHC-II+).
Mesa 3.x compatible.
"""

import mesa
import numpy as np


class CAF(mesa.Agent):
    """
    Fibroblasto asociado al cáncer (CAF).
    Subtipos basados en Öhlund et al, JEM 2017 / Elyada et al, Cancer Discov 2019:
    - myCAF: αSMA+, cerca del tumor, produce colágeno y ECM (desmoplasia)
    - iCAF: lejos del tumor, secreta IL-6, CXCL12, citoquinas inflamatorias
    - apCAF: MHC-II+, puede presentar antígenos (rol inmunomodulador)
    """

    SUBTYPES = ['myCAF', 'iCAF', 'apCAF']

    def __init__(self, model, pos=None, subtype=None):
        super().__init__(model)
        self.cell_type = 'caf'

        # Asignar subtipo (distribución realista en PDAC)
        if subtype is None:
            r = self.random.random()
            if r < 0.45:
                self.subtype = 'myCAF'    # ~45% — producen ECM
            elif r < 0.85:
                self.subtype = 'iCAF'     # ~40% — inflamatorios
            else:
                self.subtype = 'apCAF'    # ~15% — antigen presenting
        else:
            self.subtype = subtype

        self.activated = False
        self.activation_level = 0.0

    def step(self):
        # Activación por TGF-β del microambiente
        if hasattr(self.model, 'microenv') and self.pos is not None:
            me = self.model.microenv.get_local_conditions(self.pos[0], self.pos[1])
            tgfb = me.get('tgfb', 0)
            il6 = me.get('il6', 0)

            if tgfb > 0.2 or il6 > 0.2:
                self.activated = True
                self.activation_level = min(self.activation_level + 0.05, 1.0)

            # TGF-β alto → transición a myCAF (Biffi et al, Cancer Discov 2019)
            if tgfb > 0.5 and self.subtype == 'iCAF':
                if self.random.random() < 0.02:
                    self.subtype = 'myCAF'

            # IL-6 alto → transición a iCAF
            if il6 > 0.5 and self.subtype == 'myCAF':
                if self.random.random() < 0.01:
                    self.subtype = 'iCAF'

        if not self.activated:
            return

        # Comportamiento según subtipo
        if self.subtype == 'myCAF':
            self._myCAF_behavior()
        elif self.subtype == 'iCAF':
            self._iCAF_behavior()
        elif self.subtype == 'apCAF':
            self._apCAF_behavior()

        # Movimiento lento
        if self.random.random() < 0.05:
            self._move()

    def _myCAF_behavior(self):
        """myCAF: α-SMA+, produce colágeno, fibronectina → desmoplasia."""
        if self.pos is None or not hasattr(self.model, 'microenv'):
            return
        me = self.model.microenv
        x, y = self.pos

        # Depositar ECM (colágeno I/III, fibronectina)
        if 0 <= x < me.width and 0 <= y < me.height:
            me.ecm_density[x, y] = min(me.ecm_density[x, y] + 0.06, 1.0)
            me.ecm_stiffness[x, y] = min(me.ecm_stiffness[x, y] + 0.04, 1.0)
            # Secretar TGF-β (autocrino y paracrino)
            me.tgfb[x, y] = min(me.tgfb[x, y] + 0.02, 1.0)

    def _iCAF_behavior(self):
        """iCAF: secretan IL-6, CXCL12, LIF → inflamación protumoral."""
        if self.pos is None or not hasattr(self.model, 'microenv'):
            return
        me = self.model.microenv
        x, y = self.pos

        if 0 <= x < me.width and 0 <= y < me.height:
            # Secretar IL-6 (activa JAK/STAT3 en células tumorales)
            me.il6[x, y] = min(me.il6[x, y] + 0.04, 1.0)
            # Secretar CXCL12 (quimiocina, recluta Tregs y M2)
            me.cxcl12[x, y] = min(me.cxcl12[x, y] + 0.03, 1.0)
            # Menos ECM que myCAF
            me.ecm_density[x, y] = min(me.ecm_density[x, y] + 0.01, 1.0)

    def _apCAF_behavior(self):
        """apCAF: MHC-II+, modulan inmunidad (pueden activar Tregs paradójicamente)."""
        if self.pos is None or not hasattr(self.model, 'microenv'):
            return
        me = self.model.microenv
        x, y = self.pos

        if 0 <= x < me.width and 0 <= y < me.height:
            # Secretan menos citoquinas
            me.il6[x, y] = min(me.il6[x, y] + 0.01, 1.0)
            # apCAFs pueden paradójicamente activar Tregs (Elyada 2019)
            me.tgfb[x, y] = min(me.tgfb[x, y] + 0.01, 1.0)

    def _move(self):
        if self.pos is None:
            return
        nbrs = self.model.grid.get_neighborhood(
            self.pos, moore=True, include_center=False)
        empty = [p for p in nbrs if self.model.grid.is_cell_empty(p)]
        if empty:
            self.model.grid.move_agent(self, self.random.choice(empty))
