"""
Endpoints clínicos para PDAC Digital Twin.
RECIST 1.1, CA 19-9, Performance Status (ECOG), metastasis risk,
Kaplan-Meier survival estimation.
"""

import numpy as np


# ═══════════════════════════════════════
#  RECIST 1.1 (simplificado)
# ═══════════════════════════════════════
def classify_recist(baseline_count, current_count):
    """
    Clasifica respuesta tumoral según RECIST 1.1 simplificado.
    Usa conteo celular como proxy de diámetro tumoral.
    Returns: (category, label, color, emoji)
    """
    if baseline_count == 0:
        return ('N/A', 'Sin datos', '#888888', '❓')

    change_pct = ((current_count - baseline_count) / baseline_count) * 100

    if current_count <= max(baseline_count * 0.05, 2):
        return ('CR', 'Respuesta Completa', '#00c853', '🟢')
    elif change_pct <= -30:
        return ('PR', 'Respuesta Parcial', '#2979ff', '🔵')
    elif change_pct >= 20:
        return ('PD', 'Progresión', '#ff1744', '🔴')
    else:
        return ('SD', 'Enfermedad Estable', '#ffc107', '🟡')


def get_recist_history(history):
    """
    Calcula RECIST a lo largo del tiempo.
    Returns list of (step, category, change_pct).
    """
    if not history.get('cancer_alive') or len(history['cancer_alive']) < 2:
        return []

    baseline = history['cancer_alive'][0]
    results = []
    for i, count in enumerate(history['cancer_alive']):
        cat, label, color, emoji = classify_recist(baseline, count)
        change_pct = ((count - baseline) / max(baseline, 1)) * 100
        results.append({
            'step': history['step'][i],
            'category': cat,
            'label': label,
            'change_pct': change_pct,
            'color': color,
            'emoji': emoji,
        })
    return results


# ═══════════════════════════════════════
#  CA 19-9 (Biomarcador sérico)
# ═══════════════════════════════════════
CA199_NORMAL_UPPER = 37  # U/mL (límite superior normal)

def estimate_ca199(alive_count, baseline_count=30, mutation_profile=None):
    """
    Estima CA 19-9 basado en carga tumoral.
    CA 19-9 es producido por ~85% de PDAC.
    Rango clínico: normal <37 U/mL, PDAC típico 100-10,000+ U/mL.

    Modelo: CA 19-9 ∝ tumor_burden^0.8 (no lineal en tumores grandes)
    """
    # ~5-10% de pacientes son Lewis antigen-negative → no producen CA 19-9
    if mutation_profile and not mutation_profile.get('KRAS', True):
        # Sin KRAS → producción reducida
        secretion = 0.5
    else:
        secretion = 1.0

    if alive_count <= 0:
        return CA199_NORMAL_UPPER * 0.5  # Remisión

    # Modelo: baseline tumoral produce ~500 U/mL
    # Escala con tumor_burden^0.8 (sublineal para tumores grandes)
    ratio = alive_count / max(baseline_count, 1)
    ca199 = 500 * (ratio ** 0.8) * secretion

    # Ruido biológico (±15%)
    ca199 *= np.random.uniform(0.85, 1.15)

    return max(ca199, 5.0)  # Mínimo detectable


def get_ca199_trend(history, baseline_count=30, mutation_profile=None):
    """Calcula CA 19-9 a lo largo del tiempo."""
    rng = np.random.default_rng(42)  # Reproducible
    values = []
    for count in history.get('cancer_alive', []):
        ratio = count / max(baseline_count, 1)
        ca199 = 500 * (ratio ** 0.8) * rng.uniform(0.85, 1.15)
        values.append(max(ca199, 5.0))
    return values


# ═══════════════════════════════════════
#  ECOG Performance Status
# ═══════════════════════════════════════
def estimate_ecog(alive_count, baseline_count=30, treatment_toxicity=0.0):
    """
    Estima ECOG Performance Status (0-4):
    0: Asintomático
    1: Sintomático pero ambulatorio
    2: En cama <50% del día
    3: En cama >50% del día
    4: Encamado
    """
    if alive_count <= 0:
        return 0  # Sin tumor → asintomático

    ratio = alive_count / max(baseline_count, 1)

    # Score base por carga tumoral
    if ratio < 0.3:
        score = 0
    elif ratio < 0.8:
        score = 1
    elif ratio < 1.5:
        score = 2
    elif ratio < 3.0:
        score = 3
    else:
        score = 4

    # Toxicidad del tratamiento puede empeorar ECOG
    if treatment_toxicity > 0.5:
        score = min(score + 1, 4)

    return score


ECOG_DESCRIPTIONS = {
    0: "Asintomático — actividad normal",
    1: "Sintomático — ambulatorio, capaz de trabajo ligero",
    2: "Ambulatorio >50% del día — autocuidado, no trabajo",
    3: "En cama >50% del día — autocuidado limitado",
    4: "Encamado — asistencia total necesaria",
}


# ═══════════════════════════════════════
#  RIESGO DE METÁSTASIS
# ═══════════════════════════════════════
def estimate_metastasis_risk(history, mutation_profile=None):
    """
    Estima riesgo acumulativo de metástasis basado en:
    - EMT (basal-like %)
    - Señal de invasión
    - Duración del tumor
    - Angiogénesis (VEGF)

    En PDAC, ~80% tienen metástasis al diagnóstico.
    Returns: risk score [0-100%]
    """
    if not history.get('basal_like_pct') or len(history['step']) < 2:
        return 0.0

    # Componentes de riesgo
    avg_basal = np.mean(history['basal_like_pct'][-20:]) if len(history['basal_like_pct']) >= 20 \
        else np.mean(history['basal_like_pct'])

    duration_hours = history['step'][-1]
    duration_months = duration_hours / (24 * 30)

    # Tumor burden actual vs inicial
    if history['cancer_alive'][0] > 0:
        burden_ratio = history['cancer_alive'][-1] / history['cancer_alive'][0]
    else:
        burden_ratio = 0

    # Componentes
    risk = 0.0
    risk += min(avg_basal * 0.4, 40)        # EMT → metástasis (hasta 40%)
    risk += min(duration_months * 8, 30)      # Tiempo (hasta 30%)
    risk += min(max(burden_ratio - 1, 0) * 15, 20)  # Carga tumoral (hasta 20%)

    # SMAD4 loss → metástasis difusas (Iacobuzio-Donahue 2009)
    if mutation_profile and mutation_profile.get('SMAD4', False):
        risk += 10

    return np.clip(risk, 0, 95)


# ═══════════════════════════════════════
#  KAPLAN-MEIER ESTIMATION
# ═══════════════════════════════════════
def run_kaplan_meier(model_class, n_simulations=10, steps=4320,
                     protocol_name=None, mutation_profile=None,
                     model_kwargs=None, event_threshold=0.1):
    """
    Ejecuta N simulaciones con semillas distintas y estima PFS/OS.
    event_threshold: fracción del tumor inicial a la que se considera 'evento'
    Returns: dict con survival_times, median_pfs, km_curve
    """
    from drugs.treatment_protocols import get_active_drugs_for_hour

    kwargs = model_kwargs or {}
    survival_times = []
    all_curves = []

    for i in range(n_simulations):
        model = model_class(seed=i * 137, mutation_profile=mutation_profile, **kwargs)

        baseline = None
        event_step = None
        curve = []

        for step in range(steps):
            # Protocol-based dosing
            if protocol_name:
                doses = get_active_drugs_for_hour(protocol_name, step)
                model.set_drug_doses(doses)

            model.step()

            alive = model.history['cancer_alive'][-1]
            if baseline is None:
                baseline = alive

            curve.append(alive)

            # PFS event: tumor grows >20% from nadir (RECIST PD)
            nadir = min(model.history['cancer_alive'])
            if event_step is None and alive > nadir * 1.2 and step > 48:
                event_step = step

        if event_step is None:
            event_step = steps  # Censored

        survival_times.append(event_step)
        all_curves.append(curve)

    # Kaplan-Meier curve
    survival_times = sorted(survival_times)
    n = len(survival_times)
    km_times = [0]
    km_survival = [1.0]

    for i, t in enumerate(survival_times):
        km_times.append(t)
        km_survival.append(1.0 - (i + 1) / n)

    # Median PFS
    median_pfs = np.median(survival_times)

    return {
        'survival_times': survival_times,
        'median_pfs_hours': median_pfs,
        'median_pfs_days': median_pfs / 24,
        'km_times': km_times,
        'km_survival': km_survival,
        'n_simulations': n_simulations,
        'all_curves': all_curves,
    }
