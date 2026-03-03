"""
Diseñador de fármacos para nuevas dianas.
Genera moléculas hipotéticas con SMILES, analiza propiedades
y permite explorar inhibición de cualquier proteína del modelo.

Autor: Gemelo Digital PDAC 2026
"""

import numpy as np

# Intentar importar RDKit (opcional — funciona sin él con SMILES predefinidos)
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Draw, AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False


# =============================================
#  PLANTILLAS DE SMILES PARA CADA DIANA
# =============================================
# Estos SMILES son hipotéticos / inspirados en inhibidores reales
SMILES_TEMPLATES = {
    'KRAS_active': {
        'scaffold': 'c1ccc2c(c1)nc(n2)N1CCN(CC1)C(=O)c1cnc2ccccc2n1',
        'name': 'KRASi-GD01',
        'class': 'Inhibidor covalente KRAS G12D',
        'inspiration': 'MRTX1133 / Sotorasib scaffold',
        'mechanism': (
            'Se une al bolsillo Switch-II de KRAS G12D, '
            'bloqueando la interacción con efectores RAF/PI3K. '
            'Selectividad por la mutación G12D vs G12C.'
        ),
    },
    'YAP_nuclear': {
        'scaffold': 'COc1ccc(-c2nc3cc(F)ccc3o2)cc1O',
        'name': 'YAPi-GD01',
        'class': 'Disruptor de interacción YAP-TEAD',
        'inspiration': 'Verteporfin / TEADi scaffold',
        'mechanism': (
            'Bloquea la translocación nuclear de YAP '
            'e impide su unión al factor transcripcional TEAD. '
            'Reduce programas de supervivencia y resistencia a KRASi.'
        ),
    },
    'TEAD_active': {
        'scaffold': 'CC(=O)Nc1ccc(-c2ccc(F)c(NC(=O)c3ccncc3)c2)cc1',
        'name': 'TEADi-GD01',
        'class': 'Inhibidor alostérico de TEAD',
        'inspiration': 'VT-104 / K-975 scaffold',
        'mechanism': (
            'Se une al bolsillo lipofilico de TEAD, '
            'impidiendo la palmitoilación necesaria para su actividad. '
            'Bloquea la transcripción de genes pro-supervivencia.'
        ),
    },
    'HIF2A_active': {
        'scaffold': 'Oc1ccc(-c2nc(-c3ccc(F)cc3)no2)c(O)c1',
        'name': 'HIF2Ai-GD01',
        'class': 'Antagonista de HIF-2α',
        'inspiration': 'Belzutifan / PT2385 scaffold',
        'mechanism': (
            'Desestabiliza el heterodímero HIF-2α/HIF-1β, '
            'reduciendo la transcripción de genes de hipoxia. '
            'Revierte el efecto Warburg y la adicción a glutamina.'
        ),
    },
    'mTOR_active': {
        'scaffold': 'Cn1c(=O)c2c(nc(Nc3ccc(N4CCOCC4)cc3)n2C)n(C)c1=O',
        'name': 'mTORi-GD01',
        'class': 'Inhibidor dual mTORC1/mTORC2',
        'inspiration': 'Everolimus / Vistusertib scaffold',
        'mechanism': (
            'Inhibe ambos complejos mTOR, bloqueando '
            'la síntesis proteica y el crecimiento celular. '
            'Sinérgico con inhibidores de KRAS.'
        ),
    },
    'AKT_active': {
        'scaffold': 'CC(O)(c1ccc(NC(=O)c2ccncc2)cc1)c1ccc(F)cc1',
        'name': 'AKTi-GD01',
        'class': 'Inhibidor alostérico de AKT',
        'inspiration': 'MK-2206 / Ipatasertib scaffold',
        'mechanism': (
            'Se une al dominio PH de AKT, impidiendo su '
            'reclutamiento a la membrana y activación. '
            'Restaura la señal pro-apoptótica.'
        ),
    },
    'ERK_active': {
        'scaffold': 'Cc1cc(NC(=O)c2ccc(F)c(NC(=O)C3CC3)c2)no1',
        'name': 'ERKi-GD01',
        'class': 'Inhibidor de ERK1/2',
        'inspiration': 'Ulixertinib / SCH772984 scaffold',
        'mechanism': (
            'Inhibe ERK1/2, la kinasa terminal de la cascada RAS-MAPK. '
            'Útil cuando hay reactivación de ERK tras KRASi. '
            'Bloquea proliferación y supervivencia.'
        ),
    },
    'PI3K_active': {
        'scaffold': 'Cc1nc2cnc(N)nc2n1-c1ccc(S(=O)(=O)N(C)C)cc1',
        'name': 'PI3Ki-GD01',
        'class': 'Inhibidor selectivo de PI3Kα',
        'inspiration': 'Alpelisib scaffold',
        'mechanism': (
            'Inhibe selectivamente PI3Kα, cortando la señalización '
            'de supervivencia AKT/mTOR. Particularmente eficaz '
            'en células con bypass KRAS→PI3K.'
        ),
    },
    'autophagy': {
        'scaffold': 'CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12',
        'name': 'AutoI-GD01',
        'class': 'Inhibidor de autofagia (tipo cloroquina)',
        'inspiration': 'Hydroxychloroquine scaffold',
        'mechanism': (
            'Inhibe la fusión autofagosoma-lisosoma, '
            'bloqueando el reciclaje de nutrientes. '
            'En PDAC con KRAS mutado, la autofagia es esencial '
            'para la supervivencia bajo estrés metabólico.'
        ),
    },
    'macropinocytosis': {
        'scaffold': 'CC(=O)Nc1ccc(NC(=O)c2ccc(-c3ccccn3)cc2)cc1',
        'name': 'MacroI-GD01',
        'class': 'Inhibidor de macropinocitosis',
        'inspiration': 'EIPA / amiloride analogue',
        'mechanism': (
            'Bloquea la macropinocitosis mediada por KRAS oncogénico, '
            'impidiendo la captación no selectiva de proteínas extracelulares. '
            'Reduce el suministro de aminoácidos al tumor.'
        ),
    },
    'PDL1_expression': {
        'scaffold': 'CC(C)(O)c1ccc(-c2cnc3c(N)ncnc3c2)cc1',
        'name': 'PDL1i-GD01',
        'class': 'Inhibidor de molécula pequeña PD-L1',
        'inspiration': 'BMS-986189 scaffold',
        'mechanism': (
            'Molécula pequeña que induce dimerización de PD-L1, '
            'bloqueando su interacción con PD-1. '
            'Alternativa oral a anticuerpos anti-PD-L1.'
        ),
    },
}

# Propiedades de Lipinski para evaluar druglikeness
LIPINSKI_THRESHOLDS = {
    'MW': 500,      # Peso molecular < 500
    'LogP': 5,      # Lipofilicidad < 5
    'HBD': 5,       # Donadores de H-bond < 5
    'HBA': 10,      # Aceptores de H-bond < 10
}


def design_drug_for_target(target: str, inhibition_pct: float = 80) -> dict:
    """
    Diseña una molécula hipotética para la diana seleccionada.

    Args:
        target: nombre del nodo de señalización (ej: 'KRAS_active')
        inhibition_pct: porcentaje de inhibición deseado (0-100)

    Returns:
        dict con SMILES, nombre, descripción, propiedades
    """
    template = SMILES_TEMPLATES.get(target)
    if template is None:
        return {
            'success': False,
            'error': f'No hay plantilla para la diana: {target}',
            'suggestion': list(SMILES_TEMPLATES.keys()),
        }

    smiles = template['scaffold']
    result = {
        'success': True,
        'target': target,
        'inhibition_pct': inhibition_pct,
        'name': template['name'],
        'drug_class': template['class'],
        'smiles': smiles,
        'inspiration': template['inspiration'],
        'mechanism': template['mechanism'],
    }

    # --- Calcular propiedades con RDKit si está disponible ---
    if RDKIT_AVAILABLE:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            props = {
                'MW': round(Descriptors.MolWt(mol), 1),
                'LogP': round(Descriptors.MolLogP(mol), 2),
                'HBD': Descriptors.NumHDonors(mol),
                'HBA': Descriptors.NumHAcceptors(mol),
                'TPSA': round(Descriptors.TPSA(mol), 1),
                'RotBonds': Descriptors.NumRotatableBonds(mol),
                'Rings': Descriptors.RingCount(mol),
            }

            # Evaluar regla de Lipinski
            violations = 0
            if props['MW'] > LIPINSKI_THRESHOLDS['MW']:
                violations += 1
            if props['LogP'] > LIPINSKI_THRESHOLDS['LogP']:
                violations += 1
            if props['HBD'] > LIPINSKI_THRESHOLDS['HBD']:
                violations += 1
            if props['HBA'] > LIPINSKI_THRESHOLDS['HBA']:
                violations += 1

            props['lipinski_violations'] = violations
            props['druglike'] = violations <= 1

            result['properties'] = props
            result['rdkit_available'] = True
        else:
            result['properties'] = _estimate_properties(smiles)
            result['rdkit_available'] = True
            result['rdkit_parse_error'] = True
    else:
        # Estimación sin RDKit
        result['properties'] = _estimate_properties(smiles)
        result['rdkit_available'] = False

    # --- Predicción de eficacia ---
    result['predicted_ic50_nM'] = _estimate_ic50(target, inhibition_pct)
    result['selectivity_notes'] = _get_selectivity_notes(target)

    return result


def _estimate_properties(smiles: str) -> dict:
    """Estima propiedades moleculares básicas sin RDKit."""
    # Estimaciones muy aproximadas basadas en longitud SMILES
    heavy_atoms = sum(1 for c in smiles if c.isalpha() and c.isupper())
    return {
        'MW': round(heavy_atoms * 14.5 + 50, 1),  # Estimación burda
        'LogP': round(smiles.count('c') * 0.3 - smiles.count('O') * 0.5, 2),
        'HBD': smiles.count('N') + smiles.count('O') - smiles.count('n'),
        'HBA': smiles.count('N') + smiles.count('O') + smiles.count('n'),
        'lipinski_violations': 0,
        'druglike': True,
        'estimated': True,
    }


def _estimate_ic50(target: str, inhibition_pct: float) -> float:
    """Estima IC50 en nM basada en la diana y nivel de inhibición."""
    # IC50s típicos para cada clase de diana (nM)
    base_ic50 = {
        'KRAS_active': 5.0,       # MRTX1133 tiene IC50 ~0.5 nM
        'YAP_nuclear': 100.0,     # Verteporfin ~100-500 nM
        'TEAD_active': 50.0,
        'HIF2A_active': 10.0,     # Belzutifan ~10 nM
        'mTOR_active': 2.0,       # Everolimus ~2 nM
        'AKT_active': 25.0,
        'ERK_active': 5.0,
        'PI3K_active': 10.0,
        'autophagy': 500.0,       # HCQ ~1-5 µM
        'macropinocytosis': 1000.0,
        'PDL1_expression': 50.0,
    }

    base = base_ic50.get(target, 100.0)
    # Ajustar por inhibición deseada
    adjusted = base * (100.0 / max(inhibition_pct, 10.0))
    return round(adjusted, 1)


def _get_selectivity_notes(target: str) -> str:
    """Devuelve notas sobre selectividad para la diana."""
    notes = {
        'KRAS_active': (
            'Alta selectividad por G12D vs G12C vs WT. '
            'Posible toxicidad GI por inhibición de RAS en tejido normal.'
        ),
        'YAP_nuclear': (
            'Selectividad moderada. YAP es importante en regeneración tisular, '
            'posible toxicidad hepática y cutánea.'
        ),
        'TEAD_active': (
            'TEAD 1-4 tienen bolsillos similares. '
            'Un inhibidor pan-TEAD puede tener efectos off-target en corazón.'
        ),
        'HIF2A_active': (
            'Alta selectividad HIF-2α vs HIF-1α. '
            'Bien tolerado (aprobado en ccRCC como Belzutifan).'
        ),
        'mTOR_active': (
            'mTORC1 preferencial con inhibidores rapalog. '
            'Inhibidores ATP-competitivos inhiben ambos complejos.'
        ),
        'AKT_active': (
            'Tres isoformas (AKT1/2/3). '
            'Hiperglucemia es efecto adverso clase por inhibición de AKT2.'
        ),
        'ERK_active': (
            'ERK1 y ERK2 altamente homólogos. '
            'Toxicidad cutánea y ocular observada en ensayos clínicos.'
        ),
        'PI3K_active': (
            'Selectividad por isoforma α crucial para tolerar. '
            'PI3Kδ y γ son importantes en sistema inmune.'
        ),
        'autophagy': (
            'Inhibidores lisosomales poco selectivos. '
            'Toxicidad retiniana a largo plazo (cloroquina).'
        ),
        'macropinocytosis': (
            'Baja selectividad con amiloride/EIPA. '
            'Nuevos inhibidores específicos en desarrollo.'
        ),
        'PDL1_expression': (
            'Moléculas pequeñas anti-PD-L1 menos potentes que anticuerpos '
            'pero con ventaja de administración oral.'
        ),
    }
    return notes.get(target, 'Sin datos de selectividad disponibles.')


def get_combination_rationale(targets: list) -> str:
    """
    Genera un razonamiento para la combinación de dianas.

    Args:
        targets: lista de dianas (ej: ['KRAS_active', 'YAP_nuclear'])

    Returns:
        Texto con el razonamiento biológico de la combinación
    """
    rationales = {
        frozenset({'KRAS_active', 'YAP_nuclear'}): (
            '🔬 COMBINACIÓN KRASi + YAPi:\n'
            'Fundamento: YAP es el principal mecanismo de bypass tras '
            'inhibición de KRAS en PDAC. Datos preclínicos muestran que '
            'la co-inhibición previene la adquisición de resistencia y '
            'produce respuestas más duraderas.'
        ),
        frozenset({'KRAS_active', 'autophagy'}): (
            '🔬 COMBINACIÓN KRASi + Inhibidor de autofagia:\n'
            'Fundamento: Las células PDAC dependen de autofagia para '
            'sobrevivir el estrés metabólico. KRASi + HCQ mostró sinergia '
            'en modelos murinos de PDAC (ensayo HALO-301 inspiración).'
        ),
        frozenset({'KRAS_active', 'ERK_active'}): (
            '🔬 COMBINACIÓN KRASi + ERKi:\n'
            'Fundamento: Bloqueo vertical de la cascada RAS-MAPK. '
            'Previene la reactivación de ERK por feedback loops. '
            'Precaución: toxicidad aditiva cutánea y GI.'
        ),
        frozenset({'HIF2A_active', 'mTOR_active'}): (
            '🔬 COMBINACIÓN HIF2Ai + mTORi:\n'
            'Fundamento: mTOR promueve traducción de HIF-2α. '
            'Doble bloqueo del eje hipoxia-metabolismo. '
            'Eficaz en tumores hipóxicos centrales.'
        ),
    }

    key = frozenset(targets)
    if key in rationales:
        return rationales[key]

    # Generar razonamiento genérico
    target_names = ', '.join(targets)
    return (
        f'🔬 COMBINACIÓN {target_names}:\n'
        f'Combinación multi-diana que ataca {len(targets)} nodos simultáneamente. '
        f'Potencial sinérgico por bloqueo de vías compensatorias. '
        f'Se recomienda evaluar toxicidad con modelos in vitro antes de avanzar.'
    )


def get_available_targets() -> list:
    """Devuelve la lista de dianas disponibles para diseño."""
    return list(SMILES_TEMPLATES.keys())
