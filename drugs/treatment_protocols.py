"""
Protocolos de tratamiento clínicos reales para PDAC.
Cada protocolo define ciclos, fármacos, dosis y scheduling exacto.
Basado en NCCN Guidelines 2025 y ensayos publicados.
"""


PROTOCOLS = {
    # ═══════════════════════════════════════
    # PRIMERA LÍNEA
    # ═══════════════════════════════════════
    'FOLFIRINOX': {
        'name': 'FOLFIRINOX',
        'description': (
            'Oxaliplatino + Irinotecán + 5-FU + Leucovorina. '
            'Primera línea ECOG 0-1. mOS 11.1m (Conroy 2011). '
            'Más eficaz pero más tóxico que Gem/NabPac.'
        ),
        'line': '1L',
        'cycle_days': 14,
        'max_cycles': 12,
        'drugs': {
            'oxaliplatin':  {'dose': 0.85, 'days': [1], 'hours': 2},
            'irinotecan':   {'dose': 0.80, 'days': [1], 'hours': 1.5},
            '5fu':          {'dose': 0.70, 'days': [1, 2], 'hours': 46},
        },
        'reference': 'Conroy T et al. NEJM 2011;364:1817-25',
    },

    'FOLFIRINOX_mod': {
        'name': 'mFOLFIRINOX',
        'description': (
            'FOLFIRINOX modificado (sin bolus de 5-FU). '
            'Menos tóxico, similar eficacia. mOS ~10.4m.'
        ),
        'line': '1L',
        'cycle_days': 14,
        'max_cycles': 12,
        'drugs': {
            'oxaliplatin':  {'dose': 0.75, 'days': [1], 'hours': 2},
            'irinotecan':   {'dose': 0.70, 'days': [1], 'hours': 1.5},
            '5fu':          {'dose': 0.60, 'days': [1, 2], 'hours': 46},
        },
        'reference': 'Marsh R et al. Pancreas 2015;44(6):937-43',
    },

    'GemNabPac': {
        'name': 'Gemcitabina + nab-Paclitaxel',
        'description': (
            'Gemcitabina + Abraxane. Primera línea estándar. '
            'mOS 8.5m, ORR 23% (Von Hoff 2013). '
            'Mejor tolerado que FOLFIRINOX.'
        ),
        'line': '1L',
        'cycle_days': 28,
        'max_cycles': 6,
        'drugs': {
            'gemcitabine':    {'dose': 0.75, 'days': [1, 8, 15], 'hours': 0.5},
            'nab_paclitaxel': {'dose': 0.80, 'days': [1, 8, 15], 'hours': 0.5},
        },
        'reference': 'Von Hoff DD et al. NEJM 2013;369:1691-703',
    },

    'GemMono': {
        'name': 'Gemcitabina monoterapia',
        'description': (
            'Primera línea para ECOG 2 o frágiles. '
            'mOS 6.7m, ORR 7% (Burris 1997).'
        ),
        'line': '1L',
        'cycle_days': 28,
        'max_cycles': 6,
        'drugs': {
            'gemcitabine': {'dose': 0.70, 'days': [1, 8, 15], 'hours': 0.5},
        },
        'reference': 'Burris HA et al. JCO 1997;15:2403-13',
    },

    # ═══════════════════════════════════════
    # SEGUNDA LÍNEA
    # ═══════════════════════════════════════
    'NalIRI_5FU': {
        'name': 'Nal-Irinotecán + 5-FU (NAPOLI-1)',
        'description': (
            'Segunda línea tras gemcitabina. '
            'mOS 6.1m vs 4.2m (Wang-Gillam 2016).'
        ),
        'line': '2L',
        'cycle_days': 14,
        'max_cycles': 12,
        'drugs': {
            'irinotecan': {'dose': 0.65, 'days': [1], 'hours': 1.5},
            '5fu':        {'dose': 0.55, 'days': [1, 2], 'hours': 46},
        },
        'reference': 'Wang-Gillam A et al. Lancet 2016;387:545-57',
    },

    # ═══════════════════════════════════════
    # MANTENIMIENTO
    # ═══════════════════════════════════════
    'Olaparib_maint': {
        'name': 'Olaparib mantenimiento (POLO)',
        'description': (
            'Mantenimiento con PARPi tras platino en BRCA+. '
            'mPFS 7.4m vs 3.8m (Golan 2019). ~5-7% de PDAC.'
        ),
        'line': 'Maintenance',
        'cycle_days': 28,
        'max_cycles': 24,
        'drugs': {
            'olaparib': {'dose': 0.80, 'days': list(range(1, 29)), 'hours': 0},
        },
        'reference': 'Golan T et al. NEJM 2019;381:317-27',
    },

    # ═══════════════════════════════════════
    # EXPERIMENTALES / KRAS 2025-2026
    # ═══════════════════════════════════════
    'KRASi_mono': {
        'name': 'KRASi monoterapia (pan-KRAS)',
        'description': (
            'Inhibidor pan-KRAS(ON) oral diario. '
            'ORR ~20-35% en PDAC KRAS-mutado (datos fase I/II 2025).'
        ),
        'line': 'Experimental',
        'cycle_days': 28,
        'max_cycles': 12,
        'drugs': {
            'daraxonrasib': {'dose': 0.75, 'days': list(range(1, 29)), 'hours': 0},
        },
        'reference': 'Sacher A et al. NEJM 2025 (RMC-6236)',
    },

    'KRASi_combo_MEKi': {
        'name': 'KRASi + Trametinib',
        'description': (
            'Combinación pan-KRAS + MEKi para prevenir resistencia. '
            'Basado en reactivación ERK como mecanismo de escape.'
        ),
        'line': 'Experimental',
        'cycle_days': 28,
        'max_cycles': 12,
        'drugs': {
            'daraxonrasib': {'dose': 0.65, 'days': list(range(1, 29)), 'hours': 0},
            'trametinib':   {'dose': 0.50, 'days': list(range(1, 29)), 'hours': 0},
        },
        'reference': 'Preclinical rationale: Fedele C et al. Nat Med 2021',
    },

    'ICI_chemo': {
        'name': 'Anti-PD-1 + Gemcitabina',
        'description': (
            'Inmunoterapia combinada con quimio. '
            'ORR marginal en PDAC (~3-5%), excepto MSI-H.'
        ),
        'line': 'Experimental',
        'cycle_days': 21,
        'max_cycles': 8,
        'drugs': {
            'anti_pd1':    {'dose': 0.80, 'days': [1], 'hours': 0.5},
            'gemcitabine': {'dose': 0.70, 'days': [1, 8], 'hours': 0.5},
        },
        'reference': 'Weiss GJ et al. JCO 2017;35:4312',
    },

    'Sin_tratamiento': {
        'name': 'Sin tratamiento (BSC)',
        'description': 'Best supportive care. Control para comparación.',
        'line': 'Control',
        'cycle_days': 28,
        'max_cycles': 1,
        'drugs': {},
        'reference': 'N/A',
    },
}


def get_active_drugs_for_hour(protocol_name, current_hour):
    """
    Devuelve dict {drug: dose} activos en la hora actual del protocolo.
    Modela scheduling real: fármacos solo se administran en días específicos.
    """
    protocol = PROTOCOLS.get(protocol_name)
    if not protocol or not protocol['drugs']:
        return {}

    cycle_hours = protocol['cycle_days'] * 24
    max_cycles = protocol.get('max_cycles', 999)

    # ¿En qué ciclo estamos?
    cycle_num = current_hour // cycle_hours
    if cycle_num >= max_cycles:
        return {}  # Tratamiento terminado

    # Hora dentro del ciclo actual
    hour_in_cycle = current_hour % cycle_hours
    day_in_cycle = (hour_in_cycle // 24) + 1  # 1-indexed

    active = {}
    for drug_name, spec in protocol['drugs'].items():
        if day_in_cycle in spec['days']:
            # El fármaco se administra hoy
            infusion_hours = spec.get('hours', 0)
            hour_of_day = hour_in_cycle % 24
            if infusion_hours == 0:
                # Oral: activo todo el día
                active[drug_name] = spec['dose']
            elif hour_of_day < infusion_hours:
                # IV: activo solo durante infusión
                active[drug_name] = spec['dose']
            else:
                # Post-infusión: PK se encarga del decay
                active[drug_name] = spec['dose']

    return active


def get_protocol_names():
    """Devuelve lista de nombres de protocolos."""
    return list(PROTOCOLS.keys())


def get_protocol_info(name):
    """Devuelve info del protocolo."""
    return PROTOCOLS.get(name)


def get_protocol_display_names():
    """Devuelve dict {display_name: protocol_key}."""
    return {v['name']: k for k, v in PROTOCOLS.items()}
