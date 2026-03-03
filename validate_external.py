#!/usr/bin/env python3
"""
validate_external.py — Validación Externa PDAC Gemelo Digital v9

Contrasta las predicciones del modelo contra experimentos publicados
usando líneas PDAC reales (PANC-1, MiaPaCa-2, AsPC-1, BxPC-3, Capan-1).

METODOLOGÍA
──────────────────────────────────────────────────────────────────────
El modelo usa un sistema PK/PD de dos capas:

  1) El usuario pasa dose_fraction ∈ [0, 1] a set_drug_doses().
  2) TumorModel lo convierte a plasma_µM:
        plasma_µM(t) = pk_conc(t) × c_max_library × dose_fraction
  3) DrugLibrary aplica Hill n=1 con ic50_library:
        effect(plasma) = max_eff × plasma / (ic50_library + plasma)

Por tanto la conversión correcta de "concentración publicada → dose_fraction" es:
        dose_fraction = conc_µM / c_max_library

Esto garantiza que, al pico de concentración (pk_conc=1), el plasma del
modelo coincida con la concentración del experimento publicado.

La validación usa dos métricas:
  • CUANTITATIVA: viabilidad relativa (tratado/control ×100%) en el rango
    publicado ± tolerancia de 15 puntos.
  • DIRECCIONAL: comparación relativa (A < B) entre líneas celulares o
    condiciones. Más robusta ante calibración imperfecta.

La viabilidad se calcula siempre contra un control sin tratamiento con
el mismo perfil mutacional y semilla, para normalizar el crecimiento basal.

Uso:
    cd pdac_gemelo && source venv/bin/activate
    python validate_external.py
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
from simulation.tumor_model import TumorModel

# ─────────────────────────────────────────────────────────────────────────────
#  PARÁMETROS DE LA DRUG_LIBRARY (c_max en µM) — actualizados post-recalibración
#  Obtenidos directamente de drugs/drug_library.py
# ─────────────────────────────────────────────────────────────────────────────
DRUG_CMAX = {
    "gemcitabine":    60.0,   # ic50_lib=0.50 µM
    "olaparib":       20.0,   # ic50_lib=3.0 µM  (recalibrado)
    "nab_paclitaxel": 10.0,   # ic50_lib=0.10 µM
    "daraxonrasib":    5.0,   # ic50_lib=0.05 µM (= 50 nM)
    "ceralasertib":    3.0,   # ic50_lib=0.15 µM (recalibrado)
    "rsl3":            5.0,   # ic50_lib=0.08 µM
    "erastin":        15.0,   # ic50_lib=2.0 µM  (recalibrado, c_max 10→15)
    "navitoclax":      5.0,   # ic50_lib=3.0 µM  (recalibrado, c_max 2→5)
    "5fu":            20.0,   # ic50_lib=2.0 µM
    "oxaliplatin":    20.0,   # ic50_lib=5.0 µM  (recalibrado, c_max 5→20)
    "irinotecan":      4.0,   # ic50_lib=1.5 µM
}
DRUG_IC50_LIB = {
    "gemcitabine":     0.50,
    "olaparib":        3.00,   # recalibrado (revertido a 3.0, coincide con drug_library)
    "nab_paclitaxel":  0.10,
    "daraxonrasib":    0.05,
    "ceralasertib":    0.15,   # recalibrado
    "rsl3":            0.08,
    "erastin":         2.00,   # recalibrado
    "navitoclax":      3.00,   # recalibrado
    "5fu":             2.00,
    "oxaliplatin":     5.00,   # recalibrado
    "irinotecan":      1.50,
}

# ─────────────────────────────────────────────────────────────────────────────
#  IC50 PUBLICADOS (µM) POR LÍNEA CELULAR
#  (nab_ptx y daraxonrasib en µM tras conversión desde nM)
# ─────────────────────────────────────────────────────────────────────────────
PUB_IC50 = {
    # {cell_line: {drug: IC50_µM}}
    "PANC-1":   {"gemcitabine": 0.30,  "olaparib": 12.0,  "daraxonrasib": 0.045,
                 "rsl3": 0.50,  "erastin": 10.0, "navitoclax": 2.5,
                 "ceralasertib": 0.30, "nab_paclitaxel": 0.008,
                 "5fu": 3.0, "oxaliplatin": 8.0, "irinotecan": 2.0},
    "MiaPaCa-2":{"gemcitabine": 0.08,  "olaparib": 15.0,  "daraxonrasib": 0.060,
                 "rsl3": 0.30,  "erastin": 8.0,  "navitoclax": 2.0,
                 "ceralasertib": 0.20, "nab_paclitaxel": 0.005,
                 "5fu": 2.0, "oxaliplatin": 5.0, "irinotecan": 1.5},
    "AsPC-1":   {"gemcitabine": 0.50,  "olaparib": 14.0,  "daraxonrasib": 0.040,
                 "rsl3": 0.70,  "erastin": 13.0, "navitoclax": 3.0,
                 "ceralasertib": 0.40, "nab_paclitaxel": 0.012,
                 "5fu": 8.0, "oxaliplatin": 15.0, "irinotecan": 5.0},
    "BxPC-3":   {"gemcitabine": 0.15,  "olaparib": 11.0,  "daraxonrasib": 0.999,
                 "rsl3": 0.60,  "erastin": 11.0, "navitoclax": 1.8,
                 "ceralasertib": 0.50, "nab_paclitaxel": 0.010,
                 "5fu": 2.0, "oxaliplatin": 6.0, "irinotecan": 2.0},
    "Capan-1":  {"gemcitabine": 0.40,  "olaparib":  1.5,  "daraxonrasib": 0.055,
                 "rsl3": 0.40,  "erastin": 9.0,  "navitoclax": 2.2,
                 "ceralasertib": 0.25, "nab_paclitaxel": 0.009,
                 "5fu": 4.0, "oxaliplatin": 4.0, "irinotecan": 2.5},
}

# ─────────────────────────────────────────────────────────────────────────────
#  PERFILES DE LÍNEAS CELULARES
# ─────────────────────────────────────────────────────────────────────────────
CELL_LINES = {
    "PANC-1": {
        "mutation_profile": {
            "KRAS": "HET_MUT", "KRAS_variant": "G12D",
            "TP53": "HOM_LOSS", "CDKN2A": "HOM_LOSS",
            "SMAD4": "HOM_LOSS", "YAP_amp": "WT", "BRCA": "WT",
            "PTEN": "WT", "MYC": "WT", "LKB1": "WT", "MSI_status": "MSS",
        },
        "n_cancer": 40, "grid": 60, "n_caf": 12, "n_mac": 6, "n_tc": 3,
    },
    "MiaPaCa-2": {
        "mutation_profile": {
            "KRAS": "HET_MUT", "KRAS_variant": "G12C",
            "TP53": "HOM_LOSS", "CDKN2A": "HOM_LOSS",
            "SMAD4": "WT", "YAP_amp": "WT", "BRCA": "WT",
            "PTEN": "WT", "MYC": "AMP", "LKB1": "WT", "MSI_status": "MSS",
        },
        "n_cancer": 45, "grid": 60, "n_caf": 6, "n_mac": 4, "n_tc": 2,
    },
    "AsPC-1": {
        "mutation_profile": {
            "KRAS": "HET_MUT", "KRAS_variant": "G12D",
            "TP53": "HOM_LOSS", "CDKN2A": "HOM_LOSS",
            "SMAD4": "HOM_LOSS", "YAP_amp": "AMP", "BRCA": "WT",
            "PTEN": "WT", "MYC": "WT", "LKB1": "WT", "MSI_status": "MSS",
        },
        "n_cancer": 35, "grid": 55, "n_caf": 15, "n_mac": 9, "n_tc": 3,
    },
    "BxPC-3": {
        "mutation_profile": {
            "KRAS": "WT", "KRAS_variant": "WT",
            "TP53": "HET_LOSS", "CDKN2A": "HOM_LOSS",
            "SMAD4": "HOM_LOSS", "YAP_amp": "WT", "BRCA": "WT",
            "PTEN": "WT", "MYC": "WT", "LKB1": "WT", "MSI_status": "MSS",
        },
        "n_cancer": 30, "grid": 50, "n_caf": 10, "n_mac": 7, "n_tc": 5,
    },
    "Capan-1": {
        "mutation_profile": {
            "KRAS": "HET_MUT", "KRAS_variant": "G12V",
            "TP53": "WT", "CDKN2A": "HET_LOSS",
            "SMAD4": "HOM_LOSS", "YAP_amp": "WT", "BRCA": "HOM_LOSS",
            "PTEN": "WT", "MYC": "WT", "LKB1": "WT", "MSI_status": "MSS",
        },
        "n_cancer": 28, "grid": 50, "n_caf": 18, "n_mac": 5, "n_tc": 5,
    },
}

STEPS     = 72
REPEATS   = 3
SEED_BASE = 77


def dose_frac(conc_um, drug):
    """
    Conversión correcta: concentración µM → dose_fraction para el modelo.
    dose_fraction = conc_µM / c_max_library
    En el pico de concentración, plasma_µM(t_peak) ≈ c_max × dose_fraction = conc_µM.
    """
    cmax = DRUG_CMAX.get(drug, 1.0)
    return min(float(conc_um) / cmax, 1.0)


def run_exp(cl_name, doses_model, steps=STEPS, repeats=REPEATS):
    """
    Corre 'repeats' pares (control/tratado). Devuelve viabilidad media ± SD (%).
    """
    cl = CELL_LINES[cl_name]
    mp = cl["mutation_profile"]
    viabs = []
    for i in range(repeats):
        ctrl = TumorModel(
            width=cl["grid"], height=cl["grid"],
            n_cancer=cl["n_cancer"], n_caf=cl["n_caf"],
            n_macrophage=cl["n_mac"], n_tcell=cl["n_tc"],
            seed=SEED_BASE + i, mutation_profile=mp,
        )
        for _ in range(steps):
            ctrl.step()
        ctrl_n = ctrl.history["cancer_alive"][-1]

        model = TumorModel(
            width=cl["grid"], height=cl["grid"],
            n_cancer=cl["n_cancer"], n_caf=cl["n_caf"],
            n_macrophage=cl["n_mac"], n_tcell=cl["n_tc"],
            seed=SEED_BASE + i, mutation_profile=mp,
        )
        model.set_drug_doses(doses_model)
        for _ in range(steps):
            model.step()
        treated_n = model.history["cancer_alive"][-1]

        viabs.append((treated_n / max(ctrl_n, 1)) * 100)

    return float(np.mean(viabs)), float(np.std(viabs))


# ─────────────────────────────────────────────────────────────────────────────
_results = []


def check_q(test_id, desc, ref, model_v, model_sd, lit_lo, lit_hi, tol=15):
    """Test cuantitativo: model_v ∈ [lit_lo-tol, lit_hi+tol]."""
    ok = (lit_lo - tol) <= model_v <= (lit_hi + tol)
    flag = "✅ PASS" if ok else "❌ FAIL"
    print(f"  {flag} [{test_id}] {desc}")
    print(f"         Literatura: {lit_lo}–{lit_hi}% | Modelo: {model_v:.1f} ± {model_sd:.1f}%")
    print(f"         Ref: {ref}")
    _results.append({"id": test_id, "quant": True, "lit_lo": lit_lo, "lit_hi": lit_hi,
                     "model": model_v, "sd": model_sd, "pass": ok, "desc": desc})
    return ok


def check_d(test_id, desc, ref, va, la, vb, lb, expect="a<b", margin=5):
    """Test direccional: va < vb (o >) con margen."""
    ok = (va < vb - margin) if expect == "a<b" else (va > vb + margin)
    flag = "✅ PASS" if ok else "❌ FAIL"
    print(f"  {flag} [{test_id}] {desc}")
    print(f"         {la}: {va:.1f}% | {lb}: {vb:.1f}% | Δ={abs(va-vb):.1f}%")
    print(f"         Ref: {ref}")
    _results.append({"id": test_id, "quant": False, "lit_lo": None, "lit_hi": None,
                     "model": va - vb, "sd": 0, "pass": ok, "desc": desc})
    return ok


# ═════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 72)
    print("   VALIDACIÓN EXTERNA — PDAC Gemelo Digital v9")
    print("   Contraste con experimentos publicados en la literatura")
    print("=" * 72)
    print()
    print("   CONVERSIÓN DE DOSIS:  dose_fraction = conc_µM / c_max_library")
    print("   CALIBRACIÓN DRUG_LIBRARY (ic50_lib vs IC50 publicado por línea):")
    print()
    print(f"   {'Fármaco':<16} {'ic50_lib (µM)':<14} {'c_max (µM)':<12} "
          f"{'IC50 pub PANC-1 (µM)':<22} {'Factor calib.'}")
    print("   " + "─" * 72)
    for drug in DRUG_CMAX:
        lib_ic50 = DRUG_IC50_LIB[drug]
        cmax     = DRUG_CMAX[drug]
        pub_ic50 = PUB_IC50["PANC-1"].get(drug, "—")
        factor   = f"{pub_ic50 / lib_ic50:.1f}×" if isinstance(pub_ic50, float) else "—"
        flag     = "✅" if isinstance(pub_ic50, float) and abs(pub_ic50 / lib_ic50 - 1) < 5 else "⚠️"
        print(f"   {drug:<16} {lib_ic50:<14.3f} {cmax:<12.1f} "
              f"{str(pub_ic50):<22} {factor} {flag}")
    print()
    print("   ⚠️ = ic50_lib difiere >5× del IC50 publicado (recalibración recomendada)")
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE A — CALIBRACIÓN: IC50 del modelo vs publicado
    # Objetivo: ¿qué dose_fraction da ~50% viabilidad en el modelo?
    # Se usa la conversión correcta: dose_frac(ic50_pub, drug)
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE A — Calibración gemcitabina: sensibilidad diferencial por línea")
    print()

    # A1 · PANC-1 @ IC50 publicado (0.30 µM) → ~50% viabilidad
    # Bhutia 2014 Mol Pharm; Mishra 2017 Cancer Lett
    df = dose_frac(PUB_IC50["PANC-1"]["gemcitabine"], "gemcitabine")
    v_a1, sd_a1 = run_exp("PANC-1", {"gemcitabine": df})
    check_q("A1", f"Gem 0.30 µM PANC-1 @pub-IC50 (dose_frac={df:.4f}) → ~50% viabilidad",
            "Bhutia 2014 Mol Pharm; IC50 PANC-1 ~0.25–0.35 µM a 72h",
            v_a1, sd_a1, 35, 65)
    print()

    # A2 · MiaPaCa-2 @ 0.30 µM (3.75×IC50) → <40% viabilidad
    # Moore 1998 Anticancer Res
    v_a2, sd_a2 = run_exp("MiaPaCa-2",
                           {"gemcitabine": dose_frac(0.30, "gemcitabine")})
    check_q("A2", "Gem 0.30 µM MiaPaCa-2 (3.75×IC50) → <40% viabilidad",
            "Moore 1998 Anticancer Res; IC50 MiaPaCa-2 ~0.08 µM",
            v_a2, sd_a2, 5, 40)
    print()

    # A3 · AsPC-1 @ 0.30 µM (0.6×IC50) → >50% viabilidad (resistente)
    # Tan 1994 Anticancer Res; Fujita 2010
    v_a3, sd_a3 = run_exp("AsPC-1",
                           {"gemcitabine": dose_frac(0.30, "gemcitabine")})
    check_q("A3", "Gem 0.30 µM AsPC-1 (0.6×pub-IC50) → >50% viabilidad (resistente)",
            "Tan 1994 Anticancer Res; IC50 AsPC-1 ~0.5 µM (ascites, resistente)",
            v_a3, sd_a3, 50, 90)
    print()

    # A4 · BxPC-3 @ 0.15 µM (@pub-IC50) → ~50% viabilidad
    # Bai 2011 Oncol Rep
    df4 = dose_frac(PUB_IC50["BxPC-3"]["gemcitabine"], "gemcitabine")
    v_a4, sd_a4 = run_exp("BxPC-3", {"gemcitabine": df4})
    check_q("A4", f"Gem 0.15 µM BxPC-3 @pub-IC50 (dose_frac={df4:.4f}) → ~50% viabilidad",
            "Bai 2011 Oncol Rep; IC50 BxPC-3 ~0.12–0.18 µM a 72h",
            v_a4, sd_a4, 35, 65)
    print()

    # A5 · MiaPaCa-2 < AsPC-1: mayor sensibilidad a gem 0.30 µM
    # Olivares 2020 Cancers (comparativa de líneas)
    check_d("A5", "MiaPaCa-2 más sensible a gem que AsPC-1 @ 0.30 µM",
            "Olivares 2020 Cancers; IC50 MiaPaCa-2 (0.08) < AsPC-1 (0.50)",
            v_a2, "MiaPaCa-2", v_a3, "AsPC-1")
    print()

    # A6 · Capan-1 @ pub-IC50 (0.40 µM) → ~50% viabilidad
    # Bhutia 2014
    df6 = dose_frac(PUB_IC50["Capan-1"]["gemcitabine"], "gemcitabine")
    v_a6, sd_a6 = run_exp("Capan-1", {"gemcitabine": df6})
    check_q("A6", f"Gem 0.40 µM Capan-1 @pub-IC50 (dose_frac={df6:.4f}) → ~50% viabilidad",
            "Bhutia 2014 Mol Pharm; IC50 Capan-1 ~0.40 µM",
            v_a6, sd_a6, 30, 70)
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE B — OLAPARIB: letalidad sintética BRCA2 (Capan-1)
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE B — Olaparib: letalidad sintética en Capan-1 BRCA2-mutante")
    print()

    # B1 · Capan-1 @ 2×pub-IC50 (3 µM olaparib) → <45%
    # Farmer 2005 Nature; McCabe 2006 Nature
    df_b1 = dose_frac(3.0, "olaparib")
    v_b1, sd_b1 = run_exp("Capan-1", {"olaparib": df_b1})
    check_q("B1", f"Olaparib 3 µM Capan-1 (BRCA2 LOF, 2×IC50, dose_frac={df_b1:.3f}) → <45%",
            "Farmer 2005 Nature; McCabe 2006 Nature; IC50 Capan-1 ~1.5 µM",
            v_b1, sd_b1, 10, 45)
    print()

    # B2 · PANC-1 @ 3 µM olaparib (0.25×pub-IC50=12 µM) → >65%
    # Javle 2018 ASCO; BRCA2-WT no responde
    v_b2, sd_b2 = run_exp("PANC-1", {"olaparib": df_b1})  # same dose
    check_q("B2", f"Olaparib 3 µM PANC-1 (BRCA2 WT, 0.25×IC50, dose_frac={df_b1:.3f}) → >65%",
            "Javle 2018 ASCO; BRCA2-WT no responde a PARPi en monoterapia",
            v_b2, sd_b2, 65, 95)
    print()

    # B3 · Selectividad: PANC-1 − Capan-1 > 20 pp a 3 µM
    check_d("B3", "Selectividad olaparib: PANC-1 (BRCA-WT) > Capan-1 (BRCA-LOF) Δ>20%",
            "Farmer 2005 Nature; selectividad 5-15× en favor BRCA-mut",
            v_b2, "PANC-1 (BRCA-WT)", v_b1, "Capan-1 (BRCA-LOF)", expect="a>b", margin=20)
    print()

    # B4 · Escalonado de dosis Capan-1: dosis baja → efecto moderado, dosis alta → bajo
    df_low_ola  = dose_frac(1.0, "olaparib")   # baja (submaximal)
    df_high_ola = dose_frac(5.0, "olaparib")   # alta
    v_b4lo, _ = run_exp("Capan-1", {"olaparib": df_low_ola})
    v_b4hi, _ = run_exp("Capan-1", {"olaparib": df_high_ola})
    check_d("B4", "Capan-1: viabilidad disminuye con dosis creciente de olaparib",
            "Farmer 2005; McCabe 2006 — relación dosis-respuesta BRCA-def",
            v_b4lo, f"Ola 1µM ({df_low_ola:.3f})", v_b4hi, f"Ola 5µM ({df_high_ola:.3f})",
            expect="a>b", margin=5)
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE C — FERROPTOSIS: RSL3 (GPX4i) y Erastin (xCTi)
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE C — Ferroptosis: RSL3 (GPX4i) y Erastin (xCTi)")
    print()

    # C1 · RSL3 @ pub-IC50 MiaPaCa-2 (0.30 µM) → ~50%
    # Badgley 2020 Science; Dixon 2012 Cell
    df_c1 = dose_frac(PUB_IC50["MiaPaCa-2"]["rsl3"], "rsl3")
    v_c1, sd_c1 = run_exp("MiaPaCa-2", {"rsl3": df_c1})
    check_q("C1", f"RSL3 0.30 µM MiaPaCa-2 @pub-IC50 (dose_frac={df_c1:.3f}) → ~50%",
            "Badgley 2020 Science; Dixon 2012 Cell; subtipo mesenquimal hipersensible",
            v_c1, sd_c1, 20, 70)
    print()

    # C2 · Dosis-respuesta RSL3 PANC-1: IC50_lib (0.08 µM) → ~50%, 5×IC50 → <30%
    # Lim 2019 J Clin Invest
    df_c2_lo = dose_frac(DRUG_IC50_LIB["rsl3"], "rsl3")          # @ lib IC50
    df_c2_hi = dose_frac(5 * DRUG_IC50_LIB["rsl3"], "rsl3")      # @ 5×lib IC50
    v_c2lo, _ = run_exp("PANC-1", {"rsl3": df_c2_lo})
    v_c2hi, _ = run_exp("PANC-1", {"rsl3": df_c2_hi})
    check_d("C2", "RSL3 dosis-respuesta PANC-1: 5×IC50 > 1×IC50 en efecto",
            "Lim 2019 J Clin Invest; relación dosis-respuesta GPX4i",
            v_c2lo, f"@IC50_lib ({DRUG_IC50_LIB['rsl3']} µM)", v_c2hi, "@5×IC50_lib",
            expect="a>b", margin=5)
    print()

    # C3 · MiaPaCa-2 más sensible a RSL3 que AsPC-1 @ mismo dose_frac
    # Lim 2019 — KRAS G12D+SMAD4null → NRF2↑ → mayor resistencia a RSL3
    df_c3 = dose_frac(DRUG_IC50_LIB["rsl3"] * 2, "rsl3")  # 2×lib IC50
    v_c3m, _ = run_exp("MiaPaCa-2", {"rsl3": df_c3})
    v_c3a, _ = run_exp("AsPC-1",    {"rsl3": df_c3})
    check_d("C3", "MiaPaCa-2 (G12C, SMAD4-WT) más sensible a RSL3 que AsPC-1 (G12D, SMAD4-null)",
            "Lim 2019 J Clin Invest; KRAS G12D + SMAD4-null → NRF2↑ → GPX4↑ → resistencia",
            v_c3m, "MiaPaCa-2", v_c3a, "AsPC-1", margin=3)
    print()

    # C4 · Erastin @ pub-IC50 PANC-1 (10 µM) → ~50%
    # Dixon 2012 Cell; Jiang 2015 Oncogene
    df_c4 = dose_frac(PUB_IC50["PANC-1"]["erastin"], "erastin")
    v_c4, sd_c4 = run_exp("PANC-1", {"erastin": df_c4})
    check_q("C4", f"Erastin 10 µM PANC-1 @pub-IC50 (dose_frac={df_c4:.2f}) → ~50%",
            "Dixon 2012 Cell; Jiang 2015 Oncogene; IC50 PANC-1 erastin ~10 µM",
            v_c4, sd_c4, 25, 75)
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE D — PAN-KRAS (Daraxonrasib/RMC-6236): selectividad mut vs WT
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE D — Pan-KRAS (Daraxonrasib): selectividad KRAS-mut vs KRAS-WT")
    print()

    # D1 · PANC-1 (G12D) @ pub-IC50 (45 nM = 0.045 µM) → ~50%
    # Li 2023 Cancer Cell; Fell 2020 PNAS
    df_d1 = dose_frac(PUB_IC50["PANC-1"]["daraxonrasib"], "daraxonrasib")
    v_d1, sd_d1 = run_exp("PANC-1", {"daraxonrasib": df_d1})
    check_q("D1", f"Daraxonrasib 45 nM PANC-1 (G12D, @pub-IC50, dose_frac={df_d1:.4f}) → ~50%",
            "Li 2023 Cancer Cell; IC50 PANC-1 ~45 nM (G12D)",
            v_d1, sd_d1, 30, 70)
    print()

    # D2 · BxPC-3 (KRAS WT) @ mismo dose_frac → >80% (sin diana)
    # Hallin 2022 Cancer Discov
    v_d2, sd_d2 = run_exp("BxPC-3", {"daraxonrasib": df_d1})
    check_q("D2", f"Daraxonrasib 45 nM BxPC-3 (KRAS WT, dose_frac={df_d1:.4f}) → >80%",
            "Hallin 2022 Cancer Discov; KRAS-WT no tiene diana para pan-KRASi",
            v_d2, sd_d2, 80, 100)
    print()

    # D3 · MiaPaCa-2 (G12C) @ pub-IC50 (60 nM) → ~50%
    # RevMed 2023 Phase I
    df_d3 = dose_frac(PUB_IC50["MiaPaCa-2"]["daraxonrasib"], "daraxonrasib")
    v_d3, sd_d3 = run_exp("MiaPaCa-2", {"daraxonrasib": df_d3})
    check_q("D3", f"Daraxonrasib 60 nM MiaPaCa-2 (G12C, @pub-IC50, dose_frac={df_d3:.4f}) → ~50%",
            "RevMed 2023 Phase I; G12C también diana de pan-KRAS",
            v_d3, sd_d3, 30, 70)
    print()

    # D4 · PANC-1 más sensible que BxPC-3 @ mismo dose_frac
    check_d("D4", "PANC-1 (G12D) más sensible a daraxonrasib que BxPC-3 (KRAS-WT)",
            "Li 2023; Hallin 2022 — selectividad KRAS-mut/WT >100×",
            v_d1, "PANC-1 (G12D)", v_d2, "BxPC-3 (KRAS-WT)")
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE E — COMBINACIONES: sinergias publicadas
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE E — Combinaciones: sinergias validadas en la literatura")
    print()

    # E1 · Gem + Nab-ptx en PANC-1 → combo más efectivo que cada agente solo
    # Von Hoff 2013 NEJM (MPACT trial); Frese 2012 Nat Med
    df_eg = dose_frac(PUB_IC50["PANC-1"]["gemcitabine"] * 0.3, "gemcitabine")
    df_en = dose_frac(PUB_IC50["PANC-1"]["nab_paclitaxel"] * 0.3, "nab_paclitaxel")
    v_e1g, _ = run_exp("PANC-1", {"gemcitabine": df_eg})
    v_e1n, _ = run_exp("PANC-1", {"nab_paclitaxel": df_en})
    v_e1c, sd_e1c = run_exp("PANC-1", {"gemcitabine": df_eg, "nab_paclitaxel": df_en})
    ok_e1 = v_e1c < min(v_e1g, v_e1n)
    flag = "✅ PASS" if ok_e1 else "❌ FAIL"
    print(f"  {flag} [E1] Gem + Nab-ptx más efectivo que cada agente solo (PANC-1)")
    print(f"         Gem solo: {v_e1g:.1f}% | Nab solo: {v_e1n:.1f}% | Combo: {v_e1c:.1f}%")
    print(f"         Ref: Von Hoff 2013 NEJM MPACT; Frese 2012 Nat Med")
    _results.append({"id": "E1", "quant": False, "lit_lo": None, "lit_hi": None,
                     "model": v_e1c, "sd": sd_e1c, "pass": ok_e1,
                     "desc": "Gem+Nab-ptx combo > monos PANC-1"})
    print()

    # E2 · Olaparib + Ceralasertib sinergia en Capan-1 (BRCA2 LOF)
    # Kim 2020 Clin Cancer Res; O'Neil 2016 — PARPi+ATRi HRD
    df_eo = dose_frac(PUB_IC50["Capan-1"]["olaparib"] * 0.5, "olaparib")
    df_ec = dose_frac(PUB_IC50["Capan-1"]["ceralasertib"] * 0.5, "ceralasertib")
    v_e2o, _ = run_exp("Capan-1", {"olaparib": df_eo})
    v_e2c_, _ = run_exp("Capan-1", {"ceralasertib": df_ec})
    v_e2c, sd_e2c = run_exp("Capan-1", {"olaparib": df_eo, "ceralasertib": df_ec})
    ok_e2 = v_e2c < min(v_e2o, v_e2c_)
    flag = "✅ PASS" if ok_e2 else "❌ FAIL"
    print(f"  {flag} [E2] Olaparib + Ceralasertib sinergia en Capan-1 (BRCA2 LOF)")
    print(f"         Ola solo: {v_e2o:.1f}% | Cer solo: {v_e2c_:.1f}% | Combo: {v_e2c:.1f}%")
    print(f"         Ref: Kim 2020 Clin Cancer Res; O'Neil 2016 — PARPi+ATRi HRD")
    _results.append({"id": "E2", "quant": False, "lit_lo": None, "lit_hi": None,
                     "model": v_e2c, "sd": sd_e2c, "pass": ok_e2,
                     "desc": "Ola+Cer combo Capan-1"})
    print()

    # E3 · Gem baja dosis + Navitoclax (senolítico) en PANC-1
    # Ruscetti 2020 Cell — TIS priming + BCL-2i
    df_e3g = dose_frac(PUB_IC50["PANC-1"]["gemcitabine"] * 0.1, "gemcitabine")  # sub-IC50
    df_e3n = dose_frac(DRUG_IC50_LIB["navitoclax"] * 0.5, "navitoclax")          # submaximal
    v_e3g, _ = run_exp("PANC-1", {"gemcitabine": df_e3g})
    v_e3n, _ = run_exp("PANC-1", {"navitoclax": df_e3n})
    v_e3c, sd_e3c = run_exp("PANC-1", {"gemcitabine": df_e3g, "navitoclax": df_e3n})
    ok_e3 = v_e3c < min(v_e3g, v_e3n)
    flag = "✅ PASS" if ok_e3 else "❌ FAIL"
    print(f"  {flag} [E3] Gem-baja + Navitoclax (senolítico) > monos PANC-1")
    print(f"         Gem-baja solo: {v_e3g:.1f}% | Nav solo: {v_e3n:.1f}% | Combo: {v_e3c:.1f}%")
    print(f"         Ref: Ruscetti 2020 Cell; TIS priming + BCL-2i")
    _results.append({"id": "E3", "quant": False, "lit_lo": None, "lit_hi": None,
                     "model": v_e3c, "sd": sd_e3c, "pass": ok_e3,
                     "desc": "Gem-low+Nav senolytic PANC-1"})
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE F — NAVITOCLAX: calibración y diferencial entre líneas
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE F — Navitoclax: calibración y diferencial KRAS-mut vs KRAS-WT")
    print()

    # F1 · Navitoclax @ pub-IC50 PANC-1 (2.5 µM) → ~50%
    # Cang 2015 JCRCO — BCL-2/XL dependencia; ic50_lib=3.0 µM, c_max=5 µM
    df_f1 = dose_frac(PUB_IC50["PANC-1"]["navitoclax"], "navitoclax")
    v_f1, sd_f1 = run_exp("PANC-1", {"navitoclax": df_f1})
    check_q("F1", f"Navitoclax 2.5 µM PANC-1 @pub-IC50 (dose_frac={df_f1:.3f}) → ~50%",
            "Cang 2015 JCRCO; BCL-2/XL dependencia PDAC; IC50_pub PANC-1 = 2.5 µM",
            v_f1, sd_f1, 30, 70)
    print()

    # F2 · PANC-1 más sensible que BxPC-3 a navitoclax @ mismo dose_frac
    # KRAS-mut → BCL-XL dependencia elevada → más sensible a BCL-2/BCL-XLi
    v_f2, sd_f2 = run_exp("BxPC-3", {"navitoclax": df_f1})
    check_d("F2", "PANC-1 (KRAS-mut/G12D) más sensible a navitoclax que BxPC-3 (KRAS-WT)",
            "KRAS-mut → BCL-XL dependencia alta → más sensible a BCL-2/BCL-XLi",
            v_f1, "PANC-1 (G12D)", v_f2, "BxPC-3 (KRAS-WT)")
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE G — CERALASERTIB (ATRi): HRD vs HRP
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE G — Ceralasertib (ATRi): HRD (Capan-1) vs HRP (PANC-1)")
    print()

    # G1 · Capan-1 @ pub-IC50 (0.25 µM) → ~50%
    # Kim 2020 Clin Cancer Res
    df_g1 = dose_frac(PUB_IC50["Capan-1"]["ceralasertib"], "ceralasertib")
    v_g1, sd_g1 = run_exp("Capan-1", {"ceralasertib": df_g1})
    check_q("G1", f"Ceralasertib 0.25 µM Capan-1 (HRD, @pub-IC50, dose_frac={df_g1:.3f}) → ~50%",
            "Kim 2020 Clin Cancer Res; IC50 Capan-1 ceralasertib ~0.25 µM",
            v_g1, sd_g1, 30, 70)
    print()

    # G2 · Capan-1 más sensible que PANC-1 a ceralasertib @ mismo dose
    # Kim 2020; Murthy 2021 — HRD sensible a ATRi
    v_g2_panc, sd_g2 = run_exp("PANC-1", {"ceralasertib": df_g1})
    check_d("G2", "Capan-1 (BRCA2 LOF/HRD) más sensible a ceralasertib que PANC-1 (HRP)",
            "Kim 2020 Clin Cancer Res; BRCA2-def → HRD → hipersensible a ATRi",
            v_g1, "Capan-1 (HRD)", v_g2_panc, "PANC-1 (HRP)")
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE H — 5-FLUOROURACILO: sensibilidad diferencial por línea
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE H — 5-Fluorouracilo: sensibilidad diferencial por línea")
    print()

    # H1 · PANC-1 @ pub-IC50 (3 µM) → ~50% viabilidad
    # Deer 2010 Pancreas; Conroy 2011 NEJM (FOLFIRINOX)
    df_h1 = dose_frac(PUB_IC50["PANC-1"]["5fu"], "5fu")
    v_h1, sd_h1 = run_exp("PANC-1", {"5fu": df_h1})
    check_q("H1", f"5-FU 3 µM PANC-1 @pub-IC50 (dose_frac={df_h1:.3f}) → ~50% viabilidad",
            "Deer 2010 Pancreas; IC50 PANC-1 5-FU ~2-4 µM a 72h",
            v_h1, sd_h1, 35, 65)
    print()

    # H2 · AsPC-1 más resistente a 5-FU que MiaPaCa-2 @ mismo dose (3 µM)
    # Herreros-Villanueva 2012; IC50 AsPC-1 ~8 µM >> MiaPaCa-2 ~2 µM
    v_h2_mia, sd_h2_mia = run_exp("MiaPaCa-2", {"5fu": df_h1})
    v_h2_asp, sd_h2_asp = run_exp("AsPC-1",    {"5fu": df_h1})
    check_d("H2", "AsPC-1 más resistente a 5-FU que MiaPaCa-2 @ 3 µM",
            "Herreros-Villanueva 2012; IC50 AsPC-1 ~8 µM vs MiaPaCa-2 ~2 µM",
            v_h2_asp, "AsPC-1", v_h2_mia, "MiaPaCa-2", expect="a>b", margin=5)
    print()

    # H3 · FOLFIRINOX: 5-FU + Oxaliplatino + Irinotecán más efectivo que gem sola en PANC-1
    # Conroy 2011 NEJM (FOLFIRINOX vs Gem: 9.4 vs 6.8 meses OS → ~28% mejora)
    df_5fu_sub   = dose_frac(PUB_IC50["PANC-1"]["5fu"]        * 0.5, "5fu")
    df_oxali_sub = dose_frac(PUB_IC50["PANC-1"]["oxaliplatin"] * 0.5, "oxaliplatin")
    df_iri_sub   = dose_frac(PUB_IC50["PANC-1"]["irinotecan"]  * 0.5, "irinotecan")
    df_gem_sub   = dose_frac(PUB_IC50["PANC-1"]["gemcitabine"] * 0.3, "gemcitabine")
    v_h3_gem,  _ = run_exp("PANC-1", {"gemcitabine": df_gem_sub})
    v_h3_fox, sd_h3 = run_exp("PANC-1", {
        "5fu": df_5fu_sub, "oxaliplatin": df_oxali_sub, "irinotecan": df_iri_sub
    })
    check_d("H3", "FOLFIRINOX más efectivo que Gem sola en PANC-1 (sub-IC50 dosis)",
            "Conroy 2011 NEJM; FOLFIRINOX supera Gem en primera línea PDAC metastásico",
            v_h3_fox, "FOLFIRINOX", v_h3_gem, "Gem sola", expect="a<b", margin=5)
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE I — OXALIPLATINO: sensibilidad diferencial
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE I — Oxaliplatino: sensibilidad diferencial PANC-1 vs AsPC-1")
    print()

    # I1 · PANC-1 @ pub-IC50 (8 µM) → ~50%
    # Colucci 2002 JCO; Conroy 2011
    df_i1 = dose_frac(PUB_IC50["PANC-1"]["oxaliplatin"], "oxaliplatin")
    v_i1, sd_i1 = run_exp("PANC-1", {"oxaliplatin": df_i1})
    check_q("I1", f"Oxaliplatino 8 µM PANC-1 @pub-IC50 (dose_frac={df_i1:.3f}) → ~50%",
            "Colucci 2002 JCO; IC50 PANC-1 oxaliplatino ~6-10 µM",
            v_i1, sd_i1, 30, 70)
    print()

    # I2 · Capan-1 (BRCA2-LOF) más sensible que PANC-1 a oxaliplatino
    # Friboulet 2013 Clin Cancer Res; HRD → hipersensible a platinos
    v_i2_cap, sd_i2 = run_exp("Capan-1", {"oxaliplatin": df_i1})
    check_d("I2", "Capan-1 (BRCA2-LOF/HRD) más sensible a oxaliplatino que PANC-1 (BRCA-WT)",
            "Friboulet 2013 Clin Cancer Res; HRD → hipersensible a agentes platino",
            v_i2_cap, "Capan-1 (HRD)", v_i1, "PANC-1 (BRCA-WT)")
    print()

    # I3 · AsPC-1 más resistente a oxaliplatino que MiaPaCa-2 @ mismo dose
    # IC50 AsPC-1 ~15 µM >> MiaPaCa-2 ~5 µM (Yachida 2010)
    v_i3_mia, _ = run_exp("MiaPaCa-2", {"oxaliplatin": df_i1})
    v_i3_asp, _ = run_exp("AsPC-1",    {"oxaliplatin": df_i1})
    check_d("I3", "AsPC-1 más resistente a oxaliplatino que MiaPaCa-2 @ 8 µM",
            "Yachida 2010; IC50 AsPC-1 oxaliplatin ~15 µM vs MiaPaCa-2 ~5 µM",
            v_i3_asp, "AsPC-1", v_i3_mia, "MiaPaCa-2", expect="a>b", margin=3)
    print()

    # ══════════════════════════════════════════════════════════════════════
    # BLOQUE J — IRINOTECÁN: actividad e IC50 diferencial
    # ══════════════════════════════════════════════════════════════════════
    print("─" * 72)
    print("BLOQUE J — Irinotecán: sensibilidad diferencial por línea")
    print()

    # J1 · PANC-1 @ pub-IC50 (2 µM) → ~50%
    # Conroy 2011 NEJM; Maki 2005 JCO
    df_j1 = dose_frac(PUB_IC50["PANC-1"]["irinotecan"], "irinotecan")
    v_j1, sd_j1 = run_exp("PANC-1", {"irinotecan": df_j1})
    check_q("J1", f"Irinotecán 2 µM PANC-1 @pub-IC50 (dose_frac={df_j1:.3f}) → ~50%",
            "Conroy 2011 NEJM; IC50 PANC-1 irinotecán ~1.5-2.5 µM",
            v_j1, sd_j1, 30, 70)
    print()

    # J2 · MiaPaCa-2 más sensible a irinotecán que AsPC-1 @ mismo dose
    # IC50 MiaPaCa-2 ~1.5 µM vs AsPC-1 ~5 µM (Maki 2005)
    v_j2_mia, _ = run_exp("MiaPaCa-2", {"irinotecan": df_j1})
    v_j2_asp, _ = run_exp("AsPC-1",    {"irinotecan": df_j1})
    check_d("J2", "MiaPaCa-2 más sensible a irinotecán que AsPC-1 @ 2 µM",
            "Maki 2005 JCO; IC50 MiaPaCa-2 ~1.5 µM vs AsPC-1 ~5 µM",
            v_j2_mia, "MiaPaCa-2", v_j2_asp, "AsPC-1", margin=3)
    print()

    # J3 · FOLFIRINOX (5-FU+OHP+IRI) más activo que irinotecán solo en PANC-1
    # Conroy 2011; efectos aditivos de la triple combinación
    v_j3_iri, _ = run_exp("PANC-1", {"irinotecan": df_j1})
    v_j3_fox, sd_j3 = run_exp("PANC-1", {
        "5fu": dose_frac(PUB_IC50["PANC-1"]["5fu"] * 0.5, "5fu"),
        "oxaliplatin": dose_frac(PUB_IC50["PANC-1"]["oxaliplatin"] * 0.5, "oxaliplatin"),
        "irinotecan": dose_frac(PUB_IC50["PANC-1"]["irinotecan"] * 0.5, "irinotecan"),
    })
    check_d("J3", "FOLFIRINOX (sub-IC50) más efectivo que irinotecán solo en PANC-1",
            "Conroy 2011 NEJM; triple combo > cualquier monoterapia",
            v_j3_fox, "FOLFIRINOX", v_j3_iri, "Irinotecán solo", expect="a<b", margin=5)
    print()

    # ══════════════════════════════════════════════════════════════════════
    # RESUMEN FINAL + ANÁLISIS DE CALIBRACIÓN
    # ══════════════════════════════════════════════════════════════════════
    passes   = sum(1 for r in _results if r["pass"])
    total    = len(_results)
    q_pass   = sum(1 for r in _results if r["pass"] and r["quant"])
    q_total  = sum(1 for r in _results if r["quant"])
    d_pass   = sum(1 for r in _results if r["pass"] and not r["quant"])
    d_total  = sum(1 for r in _results if not r["quant"])
    pct      = passes / total * 100

    print("=" * 72)
    print(f"   RESULTADO: {passes}/{total} tests concordantes ({pct:.0f}%)")
    print(f"   • Cuantitativos:   {q_pass}/{q_total}")
    print(f"   • Direccionales:   {d_pass}/{d_total}")
    print()

    if pct >= 75:
        verdict = "✅ CONCORDANCIA ALTA — modelo reproduce bien las tendencias publicadas"
    elif pct >= 55:
        verdict = "⚡ CONCORDANCIA MODERADA — modelo parcialmente consistente"
    else:
        verdict = "❌ CONCORDANCIA BAJA — requiere recalibración de parámetros"
    print(f"   {verdict}")
    print()

    print("   TABLA RESUMEN:")
    print(f"   {'ID':<5} {'Descripción':<46} {'Lit.':<12} {'Modelo':<12} {'OK'}")
    print("   " + "─" * 72)
    for r in _results:
        flag = "✅" if r["pass"] else "❌"
        if r["quant"]:
            lit_s = f"{r['lit_lo']:.0f}–{r['lit_hi']:.0f}%"
            mod_s = f"{r['model']:.1f}±{r['sd']:.1f}%"
        else:
            lit_s = "direc."
            mod_s = f"Δ={r['model']:.1f}%"
        print(f"   {r['id']:<5} {r['desc'][:45]:<46} {lit_s:<12} {mod_s:<12} {flag}")

    print()
    print("   REFERENCIAS PRINCIPALES:")
    refs = [
        ("Farmer 2005",     "Nature 434:917",          "Olaparib letalidad sintética BRCA"),
        ("Von Hoff 2013",   "NEJM 369:1691",           "MPACT trial Gem+Nab-ptx"),
        ("Conroy 2011",     "NEJM 364:1817",           "FOLFIRINOX PDAC metastásico"),
        ("Badgley 2020",    "Science 368:85",          "Ferroptosis-cystine PDAC"),
        ("Li 2023",         "Cancer Cell 41:878",      "RMC-6236 pan-KRAS"),
        ("Ruscetti 2020",   "Cell 182:307",            "TIS + navitoclax senolytic"),
        ("Bhutia 2014",     "Mol Pharm 11:2334",       "Gem IC50 líneas PDAC"),
        ("Kim 2020",        "Clin Cancer Res 26:5405", "PARPi+ATRi sinergia HRD"),
        ("Lim 2019",        "J Clin Invest 129:3109",  "KRAS G12D NRF2/GPX4 ferroptosis"),
        ("Hallin 2022",     "Cancer Discov 12:2546",   "Selectividad KRASi KRAS-WT"),
        ("Olivares 2020",   "Cancers 12:3023",         "Gem sensibilidad diferencial"),
        ("Friboulet 2013",  "Clin Cancer Res 19:3923", "HRD sensibilidad a platinos"),
        ("Colucci 2002",    "J Clin Oncol 20:2754",    "Oxaliplatino PDAC"),
        ("Maki 2005",       "J Clin Oncol 23:4727",    "Irinotecán PDAC"),
    ]
    for a, j, t in refs:
        print(f"   • {a:<16} {j:<25} {t}")
    print()

    return passes, total


if __name__ == "__main__":
    passes, total = main()
    sys.exit(0 if passes / total >= 0.6 else 1)
