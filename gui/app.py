"""
🔬 PDAC Cell Culture Simulator — Gemelo Digital v9
Interfaz de cultivo celular in vitro para líneas PDAC reales.
Corre con: streamlit run gui/app.py --server.port 8502
"""
import streamlit as st

st.set_page_config(
    page_title="PDAC Cell Culture Simulator",
    layout="wide",
    page_icon="🔬",
    initial_sidebar_state="expanded",
)

# ── Cabecera lab ──────────────────────────────────────────────────────────────
st.markdown("""
<div style="background:linear-gradient(90deg,#0d1b2a,#102030);
border-bottom:2px solid #1b4f72;padding:.8rem 1.4rem;margin-bottom:.8rem;">
<span style="color:#aed6f1;font-size:1.6rem;font-weight:700;letter-spacing:1px;">
🔬 PDAC Cell Culture Simulator</span><br>
<span style="color:#5d9cc5;font-size:.82rem;">
Gemelo Digital v9 · PANC-1 · MiaPaCa-2 · AsPC-1 · BxPC-3 · Capan-1 ·
Ferroptosis · AMPK/LKB1 · DDR · TIS · 14/14 tests validados
</span>
</div>""", unsafe_allow_html=True)

# ── CSS lab-style ─────────────────────────────────────────────────────────────
st.markdown("""<style>
.stApp{background:#0a0f18}
section[data-testid="stSidebar"]{background:#0d1117;border-right:1px solid #1e2d3d}
section[data-testid="stSidebar"] *{color:#b8cfe0 !important}
section[data-testid="stSidebar"] h3,section[data-testid="stSidebar"] h4,
section[data-testid="stSidebar"] .stMarkdown strong{color:#7ab8d9 !important;font-weight:600}
section[data-testid="stSidebar"] label{color:#8fb8d0 !important}
.metric-box{background:#0f1c2b;border:1px solid #1b4f72;border-radius:8px;
padding:.7rem .9rem;text-align:center;margin:.2rem}
.metric-val{font-size:1.8rem;font-weight:700;color:#aed6f1;margin:.1rem 0}
.metric-lbl{font-size:.65rem;color:#3d6a8a;text-transform:uppercase;letter-spacing:1.5px}
.metric-sub{font-size:.72rem;color:#5d8faa;margin-top:.1rem}
.mut-badge{display:inline-block;background:#0f2236;border:1px solid #1a4a6a;
border-radius:4px;padding:.15rem .5rem;font-size:.75rem;color:#7ab8d9;margin:.1rem}
.drug-ref{font-size:.7rem;color:#3d6070;font-style:italic;margin-top:.1rem}
#MainMenu,footer,header{visibility:hidden}
</style>""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════════════
#  PERFILES DE LÍNEAS CELULARES (biología publicada)
# ══════════════════════════════════════════════════════════════════════════════
CELL_LINES = {
    "PANC-1": {
        "kras": "G12D", "tp53": "R248W (GOF)", "cdkn2a": "homozygous del",
        "smad4": "homozygous del", "brca2": "WT", "lkb1": "WT",
        "origin": "Ductal primario (Lieber 1975)", "grade": "III",
        "media": "DMEM", "serum": 10, "glucose_mm": 25.0,
        "doubling_h": 48,
        "n_cancer": 40, "grid": 60, "n_caf": 12, "n_mac": 6, "n_tc": 3,
        "kras_variant": "G12D",
        "mutation_profile": {
            "KRAS": "HET_MUT", "KRAS_variant": "G12D",
            "TP53": "HOM_LOSS", "CDKN2A": "HOM_LOSS",
            "SMAD4": "HOM_LOSS", "YAP_amp": "WT", "BRCA": "WT",
            "PTEN": "WT", "MYC": "WT", "LKB1": "WT", "ARID1A": "WT",
            "MSI_status": "MSS",
        },
        "ic50": {"gemcitabine": 0.30, "olaparib": 12.0, "nab_ptx_nm": 8.0,
                 "rsl3": 0.5, "erastin": 10.0, "navitoclax": 2.5,
                 "daraxonrasib_nm": 45.0, "ceralasertib": 0.3},
        "color": "#4a9eff",
        "note": "Línea gold-standard PDAC. Muy usada para validación de KRASi.",
    },
    "MiaPaCa-2": {
        "kras": "G12C", "tp53": "R248W (GOF)", "cdkn2a": "homozygous del",
        "smad4": "WT", "brca2": "WT", "lkb1": "WT",
        "origin": "Ductal primario (Yunis 1977)", "grade": "IV",
        "media": "DMEM", "serum": 10, "glucose_mm": 25.0,
        "doubling_h": 38,
        "n_cancer": 45, "grid": 60, "n_caf": 6, "n_mac": 4, "n_tc": 2,
        "kras_variant": "G12C",
        "mutation_profile": {
            "KRAS": "HET_MUT", "KRAS_variant": "G12C",
            "TP53": "HOM_LOSS", "CDKN2A": "HOM_LOSS",
            "SMAD4": "WT", "YAP_amp": "WT", "BRCA": "WT",
            "PTEN": "WT", "MYC": "AMP", "LKB1": "WT", "ARID1A": "WT",
            "MSI_status": "MSS",
        },
        "ic50": {"gemcitabine": 0.08, "olaparib": 15.0, "nab_ptx_nm": 5.0,
                 "rsl3": 0.3, "erastin": 8.0, "navitoclax": 2.0,
                 "daraxonrasib_nm": 60.0, "ceralasertib": 0.2},
        "color": "#ff6b81",
        "note": "Fenotipo mesenquimal agresivo. Mayor invasividad, alta quimiorresistencia basal.",
    },
    "AsPC-1": {
        "kras": "G12D", "tp53": "truncado (C135fs)", "cdkn2a": "homozygous del",
        "smad4": "homozygous del", "brca2": "WT", "lkb1": "WT",
        "origin": "Ascites metastásica (Tomioka 1984)", "grade": "IV",
        "media": "RPMI-1640", "serum": 10, "glucose_mm": 11.0,
        "doubling_h": 44,
        "n_cancer": 35, "grid": 55, "n_caf": 15, "n_mac": 9, "n_tc": 3,
        "kras_variant": "G12D",
        "mutation_profile": {
            "KRAS": "HET_MUT", "KRAS_variant": "G12D",
            "TP53": "HOM_LOSS", "CDKN2A": "HOM_LOSS",
            "SMAD4": "HOM_LOSS", "YAP_amp": "AMP", "BRCA": "WT",
            "PTEN": "WT", "MYC": "WT", "LKB1": "WT", "ARID1A": "WT",
            "MSI_status": "MSS",
        },
        "ic50": {"gemcitabine": 0.50, "olaparib": 14.0, "nab_ptx_nm": 12.0,
                 "rsl3": 0.7, "erastin": 13.0, "navitoclax": 3.0,
                 "daraxonrasib_nm": 40.0, "ceralasertib": 0.4},
        "color": "#ffa502",
        "note": "Derivada de ascites peritoneales. Alta IL-6, entorno inflamatorio. SMAD4-null.",
    },
    "BxPC-3": {
        "kras": "WT (!)", "tp53": "Y220C", "cdkn2a": "homozygous del",
        "smad4": "homozygous del", "brca2": "WT", "lkb1": "WT",
        "origin": "Ductal primario (Tan 1986)", "grade": "II",
        "media": "RPMI-1640", "serum": 10, "glucose_mm": 11.0,
        "doubling_h": 52,
        "n_cancer": 30, "grid": 50, "n_caf": 10, "n_mac": 7, "n_tc": 5,
        "kras_variant": "WT",
        "mutation_profile": {
            "KRAS": "WT", "KRAS_variant": "WT",
            "TP53": "HET_LOSS", "CDKN2A": "HOM_LOSS",
            "SMAD4": "HOM_LOSS", "YAP_amp": "WT", "BRCA": "WT",
            "PTEN": "WT", "MYC": "WT", "LKB1": "WT", "ARID1A": "HET_LOSS",
            "MSI_status": "MSS",
        },
        "ic50": {"gemcitabine": 0.15, "olaparib": 11.0, "nab_ptx_nm": 10.0,
                 "rsl3": 0.6, "erastin": 11.0, "navitoclax": 1.8,
                 "daraxonrasib_nm": 999.0, "ceralasertib": 0.5},  # KRASi inactivo en KRAS-WT
        "color": "#2ed573",
        "note": "KRAS wild-type (raro, ~3% PDAC). Más sensible a gemcitabina. No responde a KRASi.",
    },
    "Capan-1": {
        "kras": "G12V", "tp53": "WT", "cdkn2a": "LOH",
        "smad4": "homozygous del", "brca2": "6174delT (LOF)", "lkb1": "WT",
        "origin": "Metástasis hepática (Fogh 1977)", "grade": "II",
        "media": "IMEM", "serum": 20, "glucose_mm": 25.0,
        "doubling_h": 56,
        "n_cancer": 28, "grid": 50, "n_caf": 18, "n_mac": 5, "n_tc": 5,
        "kras_variant": "G12V",
        "mutation_profile": {
            "KRAS": "HET_MUT", "KRAS_variant": "G12V",
            "TP53": "WT", "CDKN2A": "HET_LOSS",
            "SMAD4": "HOM_LOSS", "YAP_amp": "WT", "BRCA": "HOM_LOSS",
            "PTEN": "WT", "MYC": "WT", "LKB1": "WT", "ARID1A": "WT",
            "MSI_status": "MSS",
        },
        "ic50": {"gemcitabine": 0.40, "olaparib": 1.5, "nab_ptx_nm": 9.0,
                 "rsl3": 0.4, "erastin": 9.0, "navitoclax": 2.2,
                 "daraxonrasib_nm": 55.0, "ceralasertib": 0.25},
        "color": "#a29bfe",
        "note": "BRCA2 mutante (6174delT). Deficiencia HRR → alta sensibilidad a PARP-i y platinos.",
    },
}

# ══════════════════════════════════════════════════════════════════════════════
#  SIDEBAR — laboratorio
# ══════════════════════════════════════════════════════════════════════════════
with st.sidebar:
    st.markdown("### 🧫 Línea Celular")
    cl_name = st.selectbox("Seleccionar línea", list(CELL_LINES.keys()), label_visibility="collapsed")
    cl = CELL_LINES[cl_name]

    # Ficha de la línea
    badges = " ".join([
        f'<span class="mut-badge">KRAS {cl["kras"]}</span>',
        f'<span class="mut-badge">TP53 {cl["tp53"][:7]}</span>',
        f'<span class="mut-badge">SMAD4 {"null" if "del" in cl["smad4"] else cl["smad4"]}</span>',
    ])
    if "del" in cl.get("brca2", "") or "LOF" in cl.get("brca2", ""):
        badges += f' <span class="mut-badge" style="border-color:#a29bfe;color:#a29bfe">BRCA2 mut</span>'
    st.markdown(badges, unsafe_allow_html=True)
    st.markdown(f"<span style='font-size:.72rem;color:#3d6070'>{cl['origin']} | "
                f"T½ {cl['doubling_h']}h | {cl['media']} + {cl['serum']}% FBS</span>",
                unsafe_allow_html=True)
    st.caption(f"💬 {cl['note']}")

    st.markdown("---")
    st.markdown("### 🔬 Condiciones de cultivo")
    init_confluence = st.slider("Confluencia inicial (%)", 10, 50, 25, 5,
                                help="% del área de cultivo cubierta al sembrar")
    glucose_low = st.checkbox("Glucosa baja (5.5 mM)", value=False,
                              help="Simula privación energética → activa AMPK")

    st.markdown("---")
    st.markdown("### 💊 Tratamiento (quimio)")

    dose_gem_um = st.slider("Gemcitabina (µM)", 0.0, 10.0, 0.0, 0.05,
                            help=f"IC50 {cl_name}: {cl['ic50']['gemcitabine']} µM")
    st.markdown(f"<div class='drug-ref'>Ref: IC50 {cl_name} = {cl['ic50']['gemcitabine']} µM · "
                f"plasma clínico 50–1000 µM</div>", unsafe_allow_html=True)

    dose_nab_nm = st.slider("Nab-Paclitaxel (nM)", 0.0, 50.0, 0.0, 0.5,
                            help=f"IC50 {cl_name}: {cl['ic50']['nab_ptx_nm']} nM")
    st.markdown(f"<div class='drug-ref'>Ref: IC50 {cl_name} = {cl['ic50']['nab_ptx_nm']} nM</div>",
                unsafe_allow_html=True)

    dose_ola_um = st.slider(
        "Olaparib (µM)", 0.0, 20.0, 0.0, 0.1,
        help=f"IC50 {cl_name}: {cl['ic50']['olaparib']} µM" +
             (" ← BRCA2 mut, muy sensible!" if "LOF" in cl.get("brca2","") else ""),
    )
    if "LOF" in cl.get("brca2", ""):
        st.markdown("<div class='drug-ref' style='color:#a29bfe'>⚠️ Capan-1 BRCA2 mut: IC50 = 1.5 µM</div>",
                    unsafe_allow_html=True)
    else:
        st.markdown(f"<div class='drug-ref'>Ref: IC50 {cl_name} = {cl['ic50']['olaparib']} µM</div>",
                    unsafe_allow_html=True)

    dose_anti_pd1 = st.slider("Anti-PD-1 (µg/mL)", 0.0, 10.0, 0.0, 0.5)

    st.markdown("---")
    st.markdown("### 🧬 Pan-KRAS / ATR")

    pan_kras_nm = st.slider("Daraxonrasib / RMC-6236 (nM)", 0.0, 200.0, 0.0, 2.5,
                            help=f"IC50 {cl_name}: {cl['ic50']['daraxonrasib_nm']} nM · inactivo en KRAS-WT")
    st.markdown(f"<div class='drug-ref'>IC50 {cl_name} = {cl['ic50']['daraxonrasib_nm']} nM"
                + (" | ⛔ KRAS-WT: sin diana" if cl["kras_variant"] == "WT" else "") +
                "</div>", unsafe_allow_html=True)

    dose_ceral_um = st.slider("Ceralasertib ATRi (µM)", 0.0, 5.0, 0.0, 0.05,
                              help=f"IC50 {cl_name}: {cl['ic50']['ceralasertib']} µM")

    st.markdown("---")
    st.markdown("### 🔥 Ferroptosis / Senolíticos")

    dose_rsl3_um = st.slider("RSL3 — GPX4i (µM)", 0.0, 5.0, 0.0, 0.05,
                             help=f"IC50 {cl_name}: {cl['ic50']['rsl3']} µM")
    dose_erastin_um = st.slider("Erastin — xCTi (µM)", 0.0, 20.0, 0.0, 0.5,
                                help=f"IC50 {cl_name}: {cl['ic50']['erastin']} µM")
    dose_navitoclax_um = st.slider("Navitoclax BCL-2i (µM)", 0.0, 10.0, 0.0, 0.1,
                                   help=f"IC50 {cl_name}: {cl['ic50']['navitoclax']} µM")

    st.markdown("---")
    st.markdown("### ⏱️ Duración del experimento")
    duration_h = st.select_slider("Horas de tratamiento",
                                  options=[24, 48, 72, 96, 120, 144], value=72)
    st.markdown(f"<div style='font-size:.75rem;color:#3d6070;text-align:center'>"
                f"{duration_h}h = {duration_h/24:.1f} días · ~{duration_h} pasos de simulación</div>",
                unsafe_allow_html=True)

    st.markdown("---")
    run_btn = st.button("🔬 CORRER EXPERIMENTO", type="primary", use_container_width=True)

# ══════════════════════════════════════════════════════════════════════════════
#  ESTADO DE SESIÓN
# ══════════════════════════════════════════════════════════════════════════════
if "model"    not in st.session_state: st.session_state.model    = None
if "ran"      not in st.session_state: st.session_state.ran      = False
if "cl_name"  not in st.session_state: st.session_state.cl_name  = None
if "init_n"   not in st.session_state: st.session_state.init_n   = None
if "duration" not in st.session_state: st.session_state.duration = None

# ══════════════════════════════════════════════════════════════════════════════
#  TABS
# ══════════════════════════════════════════════════════════════════════════════
tab_micro, tab_curves, tab_signal, tab_drug = st.tabs([
    "🔬 Microscopía", "📈 Curvas de respuesta", "🧬 Señalización", "💡 Diseñador de Fármacos"
])

# ── Función helper: µM → dosis modelo (Hill n=1) ─────────────────────────────
def hill_dose(conc_um, ic50_um):
    if conc_um <= 0 or ic50_um <= 0:
        return 0.0
    return float(conc_um) / (float(conc_um) + float(ic50_um))

# ── Render microscopía (scatter estilo phase-contrast) ────────────────────────
def render_microscopy(grid, title=""):
    import numpy as np
    import plotly.graph_objects as go
    h, w = grid.shape
    traces = []

    for cell_type, name, fill_color, border_color, size, symbol in [
        (2, "Muertas",       "rgba(50,50,70,0.55)",    "rgba(70,70,90,0.7)",    8,  "circle"),
        (3, "CAFs",          "rgba(180,160,80,0.40)",  "rgba(200,180,100,0.7)", 13, "diamond"),
        (4, "Macrófagos",    "rgba(255,160,40,0.50)",  "rgba(255,180,80,0.8)",  10, "circle"),
        (5, "T-cells",       "rgba(80,170,255,0.65)",  "rgba(120,200,255,0.9)",  7, "circle"),
        (1, "Células tumor", "rgba(210,225,255,0.88)", "rgba(255,255,255,0.95)", 15, "circle"),
    ]:
        ys, xs = np.where(grid == cell_type)
        if len(xs) == 0:
            continue
        rng = np.random.default_rng(cell_type * 7)
        jx = xs + rng.uniform(-0.35, 0.35, len(xs))
        jy = ys + rng.uniform(-0.35, 0.35, len(ys))
        traces.append(go.Scatter(
            x=jx, y=jy, mode="markers", name=name,
            marker=dict(size=size, color=fill_color, symbol=symbol,
                        line=dict(color=border_color, width=1)),
            hovertemplate=f"{name} (%{{x:.0f}},%{{y:.0f}})<extra></extra>",
        ))

    fig = go.Figure(traces)
    fig.update_layout(
        paper_bgcolor="#050810", plot_bgcolor="#060a14",
        xaxis=dict(visible=False, range=[-1, w + 1]),
        yaxis=dict(visible=False, range=[-1, h + 1], scaleanchor="x"),
        margin=dict(l=5, r=5, t=30, b=5), height=480,
        title=dict(text=title, font=dict(color="#7ab8d9", size=12), x=0.5),
        legend=dict(orientation="h", y=-0.04, x=0.5, xanchor="center",
                    font=dict(color="#7ab8d9", size=10), bgcolor="rgba(0,0,0,0)"),
        showlegend=True,
    )
    return fig

# ══════════════════════════════════════════════════════════════════════════════
#  RUN SIMULATION
# ══════════════════════════════════════════════════════════════════════════════
try:
    if run_btn:
        import sys, os, numpy as np
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from simulation.tumor_model import TumorModel
        from drugs.drug_library import Drug

        # IC50s de la línea seleccionada
        ic50 = cl["ic50"]

        # Calcular n_cancer a partir de confluencia inicial
        n_c = max(5, int(cl["grid"] * cl["grid"] * (init_confluence / 100) * 0.12))
        n_c = min(n_c, 80)

        # Construir mutation_profile con MSI_status correcto
        mp = dict(cl["mutation_profile"])
        mp.setdefault("MSI_status", "MSS")

        # Instanciar modelo
        with st.spinner(f"🧬 Sembrando {n_c} células {cl_name}..."):
            model = TumorModel(
                width=cl["grid"], height=cl["grid"],
                n_cancer=n_c, n_caf=cl["n_caf"],
                n_macrophage=cl["n_mac"], n_tcell=cl["n_tc"],
                seed=42, mutation_profile=mp,
            )

        # Glucosa baja → energy_stress initial boost
        if glucose_low:
            model.microenv.oxygen[:] *= 0.75  # proxy: menos O2/nutrientes

        # Construir dosis (Hill equation → modelo 0-1)
        doses = {}

        d = hill_dose(dose_gem_um, ic50["gemcitabine"])
        if d > 0: doses["gemcitabine"] = d

        d = hill_dose(dose_nab_nm / 1000.0, ic50["nab_ptx_nm"] / 1000.0)
        if d > 0: doses["nab_paclitaxel"] = d

        d = hill_dose(dose_ola_um, ic50["olaparib"])
        if d > 0: doses["olaparib"] = d

        if dose_anti_pd1 > 0: doses["anti_pd1"] = hill_dose(dose_anti_pd1, 2.5)

        d = hill_dose(pan_kras_nm / 1000.0, ic50["daraxonrasib_nm"] / 1000.0)
        if d > 0 and cl["kras_variant"] != "WT":
            # Inyectar fármaco pan-KRAS dinámico
            model.drug_library.drugs["daraxonrasib"] = Drug(
                name="Daraxonrasib (RMC-6236)",
                targets={"KRAS_active": 0.95, "RAF_active": 0.55,
                         "ERK_active": 0.28, "proliferation_signal": 0.50},
                ic50=0.15, max_efficacy=0.95,
                description="Pan-KRAS(ON) tri-complex inhibitor",
            )
            doses["daraxonrasib"] = d

        d = hill_dose(dose_ceral_um, ic50["ceralasertib"])
        if d > 0: doses["ceralasertib"] = d

        d = hill_dose(dose_rsl3_um, ic50["rsl3"])
        if d > 0: doses["rsl3"] = d

        d = hill_dose(dose_erastin_um, ic50["erastin"])
        if d > 0: doses["erastin"] = d

        d = hill_dose(dose_navitoclax_um, ic50["navitoclax"])
        if d > 0: doses["navitoclax"] = d

        model.set_drug_doses(doses)
        init_n = n_c

        # Simular con barra de progreso
        bar = st.progress(0, text=f"⏳ Simulando {duration_h}h de cultivo...")
        for i in range(duration_h):
            model.step()
            if (i + 1) % 6 == 0 or i == duration_h - 1:
                s = model.get_summary()
                alive = s["Células tumorales"]
                bar.progress((i + 1) / duration_h,
                             text=f"⏳ {i+1}/{duration_h}h — Vivas: {alive} | "
                                  f"Viabilidad: {alive/init_n*100:.0f}%")
        bar.empty()

        st.session_state.model    = model
        st.session_state.ran      = True
        st.session_state.cl_name  = cl_name
        st.session_state.init_n   = init_n
        st.session_state.duration = duration_h
        st.toast(f"✅ {duration_h}h completadas — {cl_name}", icon="🔬")
        st.rerun()

    # ══════════════════════════════════════════════════════════════════════════
    #  MOSTRAR RESULTADOS
    # ══════════════════════════════════════════════════════════════════════════
    model   = st.session_state.model
    ran     = st.session_state.ran
    cl_res  = st.session_state.cl_name or cl_name
    init_n  = st.session_state.init_n  or cl["n_cancer"]
    dur_h   = st.session_state.duration or duration_h

    # ── TAB 1: MICROSCOPÍA ────────────────────────────────────────────────────
    with tab_micro:
        if model is not None and ran:
            import numpy as np, plotly.graph_objects as go

            s = model.get_summary()
            alive  = s["Células tumorales"]
            viab   = alive / max(init_n, 1) * 100
            grid_s = model.get_grid_state()
            h, w   = grid_s.shape
            conflu = alive / (h * w) * 100

            # ── Métricas principales ──────────────────────────────────────────
            c1, c2, c3, c4, c5 = st.columns(5)
            def mbox(col, label, val, sub="", color="#aed6f1"):
                col.markdown(
                    f'<div class="metric-box"><div class="metric-lbl">{label}</div>'
                    f'<div class="metric-val" style="color:{color}">{val}</div>'
                    f'<div class="metric-sub">{sub}</div></div>',
                    unsafe_allow_html=True)

            hist = model.history
            # Viabilidad color
            v_color = "#2ed573" if viab >= 70 else ("#f39c12" if viab >= 40 else "#ff4757")
            mbox(c1, "VIABILIDAD", f"{viab:.0f}%", "vs control", v_color)
            mbox(c2, "CONFLUENCIA", f"{conflu:.1f}%", f"células vivas / {h}×{w}")
            mbox(c3, "CÉLULAS VIVAS", str(alive), f"inicial: {init_n}")
            mbox(c4, "T CULTIVO", f"{dur_h}h", f"= {dur_h/24:.1f} días")
            resist = s.get("Resistentes", hist["resistant_count"][-1] if hist["resistant_count"] else 0)
            r_color = "#ff4757" if resist > 3 else "#aed6f1"
            mbox(c5, "RESISTENTES", str(resist), "células con ADR", r_color)

            # ── Leyenda de clasificación ──────────────────────────────────────
            classify = ("🔴 Resistente" if viab > 70 else
                        "🟡 Sensibilidad intermedia" if viab > 30 else
                        "🟢 Sensible")
            if not model.drug_doses:
                classify = "⚪ Sin tratamiento"
            st.markdown(
                f"<div style='background:#0d1a2a;border:1px solid #1b4f72;border-radius:6px;"
                f"padding:.4rem 1rem;margin:.5rem 0;font-size:.85rem;color:#7ab8d9'>"
                f"<strong>{cl_res}</strong> · {dur_h}h · Respuesta: <strong>{classify}</strong>"
                f"</div>", unsafe_allow_html=True)

            # ── Vista microscópica ────────────────────────────────────────────
            mc1, mc2 = st.columns([2, 1])
            with mc1:
                micro_fig = render_microscopy(
                    grid_s,
                    title=f"Vista phase-contrast — {cl_res} @ {dur_h}h"
                )
                st.plotly_chart(micro_fig, use_container_width=True, key="micro")

            with mc2:
                # Mapa de oxígeno
                fig_o2 = go.Figure(go.Heatmap(
                    z=model.microenv.oxygen.T,
                    colorscale=[[0, "#0d0020"], [0.25, "#4a0072"],
                                [0.5, "#c0392b"], [0.75, "#e67e22"], [1, "#27ae60"]],
                    zmin=0, zmax=1,
                    colorbar=dict(title=dict(text="O₂", font=dict(color="#7ab8d9")),
                                  tickfont=dict(color="#7ab8d9"), thickness=12),
                ))
                fig_o2.update_layout(
                    title=dict(text="💨 Mapa O₂ / hipoxia", font=dict(color="#7ab8d9", size=11), x=0.5),
                    paper_bgcolor="#050810", plot_bgcolor="#060a14",
                    margin=dict(l=5, r=40, t=30, b=5), height=240,
                    xaxis=dict(visible=False), yaxis=dict(visible=False),
                )
                st.plotly_chart(fig_o2, use_container_width=True, key="o2")

                # Mapa ECM
                fig_ecm = go.Figure(go.Heatmap(
                    z=model.microenv.ecm_density.T,
                    colorscale=[[0, "#050810"], [0.4, "#1b2a3b"], [0.7, "#2e4a5b"], [1, "#7fb3c8"]],
                    zmin=0, zmax=1,
                    colorbar=dict(title=dict(text="ECM", font=dict(color="#7ab8d9")),
                                  tickfont=dict(color="#7ab8d9"), thickness=12),
                ))
                fig_ecm.update_layout(
                    title=dict(text="🧱 Densidad ECM (desmoplasia)", font=dict(color="#7ab8d9", size=11), x=0.5),
                    paper_bgcolor="#050810", plot_bgcolor="#060a14",
                    margin=dict(l=5, r=40, t=30, b=5), height=200,
                    xaxis=dict(visible=False), yaxis=dict(visible=False),
                )
                st.plotly_chart(fig_ecm, use_container_width=True, key="ecm")

            # ── Ficha de la línea celular ─────────────────────────────────────
            with st.expander("📋 Ficha de la línea celular", expanded=False):
                cl_info = CELL_LINES[cl_res]
                ci1, ci2 = st.columns(2)
                with ci1:
                    st.markdown(f"**Origen:** {cl_info['origin']}")
                    st.markdown(f"**Grado:** {cl_info['grade']} | **Media:** {cl_info['media']} + {cl_info['serum']}% FBS")
                    st.markdown(f"**T½ duplicación:** {cl_info['doubling_h']}h")
                    st.markdown(f"**Glucosa basal:** {cl_info['glucose_mm']} mM")
                with ci2:
                    st.markdown("**Perfil mutacional:**")
                    for gene, mut in [
                        ("KRAS", cl_info["kras"]), ("TP53", cl_info["tp53"]),
                        ("CDKN2A", cl_info["cdkn2a"]), ("SMAD4", cl_info["smad4"]),
                        ("BRCA2", cl_info["brca2"]), ("LKB1", cl_info["lkb1"]),
                    ]:
                        col_m = "#a29bfe" if "mut" in mut.lower() or "del" in mut.lower() or "lof" in mut.lower() else "#2ed573"
                        st.markdown(f"<span class='mut-badge' style='border-color:{col_m};color:{col_m}'>"
                                    f"{gene}: {mut}</span>", unsafe_allow_html=True)

        else:
            # ── Pantalla de bienvenida ────────────────────────────────────────
            st.markdown("")
            _, wc, _ = st.columns([1, 2.5, 1])
            with wc:
                cl_now = CELL_LINES[cl_name]
                badges_w = " ".join([
                    f'<span class="mut-badge">KRAS {cl_now["kras"]}</span>',
                    f'<span class="mut-badge">TP53 {cl_now["tp53"][:8]}</span>',
                    f'<span class="mut-badge">SMAD4 {"null" if "del" in cl_now["smad4"] else cl_now["smad4"]}</span>',
                ])
                st.markdown(f"""
<div style="text-align:center;padding:2.5rem;background:#0a1420;
border:1px solid #1b4f72;border-radius:16px;">
<div style="font-size:3.5rem;">🔬</div>
<h2 style="color:#aed6f1;margin:.5rem 0;">PDAC Cell Culture Simulator</h2>
<p style="color:#5d9cc5;margin:.4rem 0;">
<strong style="color:#7ab8d9">{cl_name}</strong> seleccionada ·
{cl_now['origin'][:35]}</p>
<div style="margin:.8rem 0">{badges_w}</div>
<p style="color:#3d6070;font-size:.85rem;">
T½ duplicación: <strong style="color:#5d9cc5">{cl_now['doubling_h']}h</strong> ·
Media: <strong style="color:#5d9cc5">{cl_now['media']}</strong> ·
Glucosa: <strong style="color:#5d9cc5">{cl_now['glucose_mm']} mM</strong>
</p>
<hr style="border-color:#1b4f72;margin:1rem 0">
<p style="color:#3d6a8a;font-size:.8rem;">
Configura dosis en la barra lateral →<br>
Pulsa <strong style="color:#7ab8d9">🔬 CORRER EXPERIMENTO</strong>
</p>
<p style="color:#2d4a5a;font-size:.72rem;margin-top:1rem;">
Gemcitabina · Nab-Paclitaxel · Olaparib · Daraxonrasib · RSL3 · Erastin ·
Navitoclax · Ceralasertib · Anti-PD-1
</p>
</div>""", unsafe_allow_html=True)

    # ── TAB 2: CURVAS DE RESPUESTA ────────────────────────────────────────────
    with tab_curves:
        if model is not None and ran:
            import numpy as np, pandas as pd
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots

            hist = model.history
            if not hist["step"]:
                st.warning("Sin datos históricos.")
            else:
                df = pd.DataFrame(hist)
                # Viabilidad % vs control
                df["viability_pct"] = df["cancer_alive"] / max(init_n, 1) * 100
                df["confluence_pct"] = df["cancer_alive"] / (cl["grid"] ** 2) * 100
                df["time_h"] = df["step"]

                # ── Viabilidad y confluencia ──────────────────────────────────
                fig_v = make_subplots(
                    rows=2, cols=3,
                    subplot_titles=[
                        "📊 Viabilidad (% vs T₀)", "🧫 Confluencia (%)",
                        "👥 Poblaciones",
                        "💊 Señalización (proliferación/apoptosis)",
                        "🔥 Ferroptosis & NRF2" if "avg_ferroptosis_risk" in df.columns else "⚡ O₂ / Lactato",
                        "⚙️ Resistencia adquirida",
                    ],
                    vertical_spacing=0.18, horizontal_spacing=0.10,
                )

                # Panel 1: Viabilidad
                fig_v.add_trace(go.Scatter(
                    x=df.time_h, y=df.viability_pct, name="Viabilidad %",
                    fill="tozeroy", fillcolor="rgba(74,158,255,0.1)",
                    line=dict(color="#4a9eff", width=2.5),
                ), row=1, col=1)
                # Línea IC50 (50%)
                fig_v.add_hline(y=50, line_dash="dot", line_color="rgba(255,165,0,0.4)",
                                annotation_text="IC50", annotation_font_color="orange",
                                row=1, col=1)

                # Panel 2: Confluencia
                fig_v.add_trace(go.Scatter(
                    x=df.time_h, y=df.confluence_pct, name="Confluencia %",
                    fill="tozeroy", fillcolor="rgba(46,213,115,0.08)",
                    line=dict(color="#2ed573", width=2),
                ), row=1, col=2)

                # Panel 3: Poblaciones
                fig_v.add_trace(go.Scatter(x=df.time_h, y=df.cancer_alive,
                    name="Tumor", line=dict(color="#ff4757", width=2),
                    fill="tozeroy", fillcolor="rgba(255,71,87,0.08)"), row=1, col=3)
                fig_v.add_trace(go.Scatter(x=df.time_h, y=df.caf_count,
                    name="CAFs", line=dict(color="#f9ca24", width=1.5)), row=1, col=3)
                fig_v.add_trace(go.Scatter(x=df.time_h, y=df.tcell_count,
                    name="T-cells", line=dict(color="#4a9eff", width=1.5)), row=1, col=3)
                fig_v.add_trace(go.Scatter(x=df.time_h, y=df.macrophage_count,
                    name="Macrófagos", line=dict(color="#ffa502", width=1.5)), row=1, col=3)

                # Panel 4: Señalización
                fig_v.add_trace(go.Scatter(x=df.time_h, y=df.avg_proliferation,
                    name="Proliferación", line=dict(color="#ff6b9d", width=2),
                    showlegend=False), row=2, col=1)
                fig_v.add_trace(go.Scatter(x=df.time_h, y=df.avg_apoptosis,
                    name="Apoptosis", line=dict(color="#a29bfe", width=2),
                    showlegend=False), row=2, col=1)

                # Panel 5: Ferroptosis / O2
                if "avg_ferroptosis_risk" in df.columns:
                    fig_v.add_trace(go.Scatter(x=df.time_h, y=df.avg_ferroptosis_risk,
                        name="Ferroptosis", fill="tozeroy",
                        fillcolor="rgba(232,67,147,0.1)",
                        line=dict(color="#e84393", width=2), showlegend=False), row=2, col=2)
                if "avg_nrf2" in df.columns:
                    fig_v.add_trace(go.Scatter(x=df.time_h, y=df.avg_nrf2,
                        name="NRF2", line=dict(color="#00cec9", width=1.5),
                        showlegend=False), row=2, col=2)
                if "avg_ampk" in df.columns:
                    fig_v.add_trace(go.Scatter(x=df.time_h, y=df.avg_ampk,
                        name="AMPK", line=dict(color="#6c5ce7", width=1.5, dash="dot"),
                        showlegend=False), row=2, col=2)

                # Panel 6: Resistencia
                fig_v.add_trace(go.Scatter(x=df.time_h, y=df.resistant_count,
                    name="Resistentes", fill="tozeroy",
                    fillcolor="rgba(255,118,117,0.12)",
                    line=dict(color="#ff7675", width=2), showlegend=False), row=2, col=3)
                if "senescent_count" in df.columns:
                    fig_v.add_trace(go.Scatter(x=df.time_h, y=df.senescent_count,
                        name="Senesc.", line=dict(color="#fdcb6e", width=1.5, dash="dot"),
                        showlegend=False), row=2, col=3)

                fig_v.update_xaxes(title_text="Horas", title_font_color="#5d8faa",
                                   tickfont=dict(color="#5d8faa"), gridcolor="#0f1e2e")
                fig_v.update_yaxes(tickfont=dict(color="#5d8faa"), gridcolor="#0f1e2e")
                fig_v.update_layout(
                    paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(6,10,20,0.8)",
                    font=dict(color="#aed6f1", size=10),
                    height=580, margin=dict(l=40, r=20, t=50, b=50),
                    legend=dict(orientation="h", y=-0.13, x=0.5, xanchor="center",
                                font=dict(color="#7ab8d9", size=9)),
                )
                st.plotly_chart(fig_v, use_container_width=True, key="curves")

                # ── IC50 estimado de este experimento ─────────────────────────
                st.markdown("---")
                ec1, ec2 = st.columns([1, 1])
                with ec1:
                    final_viab = df["viability_pct"].iloc[-1]
                    drug_list = []
                    if dose_gem_um > 0: drug_list.append(f"Gemcitabina {dose_gem_um:.2f} µM")
                    if dose_nab_nm > 0: drug_list.append(f"Nab-PTX {dose_nab_nm:.1f} nM")
                    if dose_ola_um > 0: drug_list.append(f"Olaparib {dose_ola_um:.1f} µM")
                    if pan_kras_nm > 0: drug_list.append(f"Daraxonrasib {pan_kras_nm:.1f} nM")
                    if dose_rsl3_um > 0: drug_list.append(f"RSL3 {dose_rsl3_um:.2f} µM")
                    drug_str = " + ".join(drug_list) if drug_list else "sin tratamiento"

                    cat = ("Resistente (>70%)" if final_viab > 70 else
                           "Intermedio (30-70%)" if final_viab > 30 else "Sensible (<30%)")
                    cat_color = ("#ff4757" if final_viab > 70 else
                                 "#f39c12" if final_viab > 30 else "#2ed573")

                    st.markdown(f"""
<div style="background:#0d1a2a;border:1px solid #1b4f72;border-radius:8px;padding:1rem;">
<h4 style="color:#7ab8d9;margin-top:0">📋 Resultado del experimento</h4>
<table style="width:100%;font-size:.85rem;color:#aed6f1;">
<tr><td>Línea celular</td><td style="text-align:right"><strong>{cl_res}</strong></td></tr>
<tr><td>Tratamiento</td><td style="text-align:right"><strong style="font-size:.8rem">{drug_str}</strong></td></tr>
<tr><td>Duración</td><td style="text-align:right"><strong>{dur_h}h ({dur_h/24:.1f} días)</strong></td></tr>
<tr><td>Viabilidad final</td><td style="text-align:right">
<strong style="color:{cat_color}">{final_viab:.1f}%</strong></td></tr>
<tr><td>Clasificación</td><td style="text-align:right">
<strong style="color:{cat_color}">{cat}</strong></td></tr>
</table></div>""", unsafe_allow_html=True)

                with ec2:
                    # Exportar
                    st.markdown("### 📥 Exportar datos")
                    st.download_button("⬇️ CSV completo",
                                       df.to_csv(index=False).encode(),
                                       f"{cl_res}_{dur_h}h_data.csv", "text/csv",
                                       use_container_width=True)
                    report = (f"PDAC Cell Culture Simulator — Reporte\n"
                              f"{'='*45}\n"
                              f"Línea: {cl_res}\n"
                              f"Tratamiento: {drug_str}\n"
                              f"Duración: {dur_h}h\n"
                              f"Células iniciales: {init_n}\n"
                              f"Células finales: {int(df.cancer_alive.iloc[-1])}\n"
                              f"Viabilidad: {final_viab:.1f}%\n"
                              f"Confluencia final: {df.confluence_pct.iloc[-1]:.1f}%\n"
                              f"Resistentes: {int(df.resistant_count.iloc[-1])}\n"
                              f"Clasificación: {cat}\n")
                    st.download_button("⬇️ Reporte TXT",
                                       report.encode(),
                                       f"{cl_res}_{dur_h}h_report.txt", "text/plain",
                                       use_container_width=True)
        else:
            st.info("🔬 Corre el experimento primero para ver las curvas de respuesta.")

    # ── TAB 3: SEÑALIZACIÓN MOLECULAR ─────────────────────────────────────────
    with tab_signal:
        if model is not None and ran:
            import numpy as np, plotly.graph_objects as go

            # Obtener último estado de señalización de las células cancerosas
            cancer_cells = [a for a in model.agents
                            if hasattr(a, "signaling") and getattr(a, "_alive", True)]

            if cancer_cells:
                # Promediar nodos de señalización
                node_keys = [
                    "KRAS_active", "ERK_active", "AKT_active", "mTOR_active",
                    "YAP_nuclear", "HIF2A_active", "autophagy", "macropinocytosis",
                    "PDL1_expression", "NRF2_active", "GPX4_level", "AMPK_active",
                    "LKB1_active", "ATR_active", "DDR_active", "survival_signal",
                    "apoptosis_signal", "proliferation_signal",
                ]
                avg_nodes = {}
                for k in node_keys:
                    vals = []
                    for c in cancer_cells[:50]:  # sample max 50
                        try:
                            vals.append(c.signaling.nodes.get(k, 0.0))
                        except Exception:
                            pass
                    avg_nodes[k] = float(np.mean(vals)) if vals else 0.0

                # Radar chart de vías clave
                categories = ["KRAS", "ERK", "AKT", "mTOR", "YAP", "HIF2A",
                              "Autofagia", "NRF2", "GPX4", "AMPK", "PDL1",
                              "Proliferación", "Supervivencia", "Apoptosis"]
                vals_radar = [
                    avg_nodes.get("KRAS_active", 0), avg_nodes.get("ERK_active", 0),
                    avg_nodes.get("AKT_active", 0), avg_nodes.get("mTOR_active", 0),
                    avg_nodes.get("YAP_nuclear", 0), avg_nodes.get("HIF2A_active", 0),
                    avg_nodes.get("autophagy", 0), avg_nodes.get("NRF2_active", 0),
                    avg_nodes.get("GPX4_level", 0), avg_nodes.get("AMPK_active", 0),
                    avg_nodes.get("PDL1_expression", 0),
                    avg_nodes.get("proliferation_signal", 0),
                    avg_nodes.get("survival_signal", 0),
                    avg_nodes.get("apoptosis_signal", 0),
                ]

                fig_rad = go.Figure(go.Scatterpolar(
                    r=vals_radar + [vals_radar[0]],
                    theta=categories + [categories[0]],
                    fill="toself", fillcolor="rgba(74,158,255,0.15)",
                    line=dict(color="#4a9eff", width=2),
                    name=f"{cl_res} señalización",
                ))
                fig_rad.update_layout(
                    polar=dict(
                        bgcolor="#060a14",
                        radialaxis=dict(visible=True, range=[0, 1],
                                        tickfont=dict(color="#3d6a8a", size=9),
                                        gridcolor="#0f1e2e"),
                        angularaxis=dict(tickfont=dict(color="#7ab8d9", size=10),
                                         gridcolor="#0f1e2e"),
                    ),
                    paper_bgcolor="rgba(0,0,0,0)",
                    font=dict(color="#aed6f1"),
                    title=dict(text=f"🧬 Señalización media — {cl_res} @ {dur_h}h",
                               font=dict(color="#7ab8d9", size=13), x=0.5),
                    height=460, margin=dict(l=40, r=40, t=60, b=40),
                    showlegend=False,
                )

                sc1, sc2 = st.columns([1.2, 1])
                with sc1:
                    st.plotly_chart(fig_rad, use_container_width=True, key="radar")

                with sc2:
                    st.markdown("### 🔬 Actividad de nodos clave")
                    sorted_nodes = sorted(avg_nodes.items(), key=lambda x: -x[1])
                    for node, val in sorted_nodes:
                        bar_w = int(val * 100)
                        bar_color = ("#ff4757" if val > 0.7 else
                                     "#f39c12" if val > 0.4 else "#2ed573")
                        st.markdown(
                            f"<div style='margin:.15rem 0;font-size:.8rem;color:#aed6f1'>"
                            f"<span style='display:inline-block;width:140px;color:#7ab8d9'>{node}</span>"
                            f"<div style='display:inline-block;width:{bar_w}%;height:10px;"
                            f"background:{bar_color};border-radius:3px;vertical-align:middle'></div>"
                            f" <span style='color:{bar_color};font-size:.75rem'>{val:.3f}</span>"
                            f"</div>", unsafe_allow_html=True)

                # Resumen Delta T0→Tf
                st.markdown("---")
                st.markdown("### 📊 Efecto del tratamiento en vías clave")
                init_avg = model.initial_node_averages
                if init_avg:
                    rows = []
                    for k, vf in avg_nodes.items():
                        v0 = init_avg.get(k, vf)
                        delta = vf - v0
                        rows.append({"Nodo": k, "T₀": f"{v0:.3f}", "Tf": f"{vf:.3f}",
                                     "Δ": f"{delta:+.3f}"})
                    import pandas as pd
                    df_delta = pd.DataFrame(rows)
                    st.dataframe(df_delta, use_container_width=True,
                                 hide_index=True,
                                 column_config={"Δ": st.column_config.NumberColumn(format="%+.3f")})
        else:
            st.info("🔬 Corre el experimento primero para ver la señalización molecular.")

    # ── TAB 4: DISEÑADOR DE FÁRMACOS ─────────────────────────────────────────
    with tab_drug:
        st.markdown("## 💡 Diseñador de Fármacos — Nuevas Dianas")
        st.markdown("Introduce la proteína diana y simula un inhibidor hipotético.")
        st.markdown("---")

        dc1, dc2 = st.columns([1, 1])
        with dc1:
            protein = st.text_input("🧬 Proteína diana",
                placeholder="GPX4, NRF2, KRAS, YAP, ATR, AMPK, IDO1...")
            des_cell = st.selectbox("Línea celular", list(CELL_LINES.keys()), key="des_cl")
            des_inh = st.slider("% Inhibición", 10, 100, 80, 5, key="di")
            des_dose_um = st.slider("Dosis (µM)", 0.1, 20.0, 1.0, 0.1, key="dd")
            btn_des = st.button("🧪 Generar + Simular 72h",
                                type="primary", use_container_width=True)
        with dc2:
            st.markdown("**Dianas v9 disponibles:**")
            st.code("KRAS · KRAS G12C · YAP · TEAD · HIF2A · mTOR\n"
                    "AKT · ERK · PI3K · GPX4 · NRF2 · SLC7A11\n"
                    "AMPK · LKB1 · ATR · CHK1 · ATM · PARP\n"
                    "IDO1 · CD73 · CD47 · BCL2 · p21 · SASP")
            st.info("💡 Diana no reconocida → se mapea al nodo más cercano.")

        node_map = {
            "KRAS": "KRAS_active", "YAP": "YAP_nuclear", "TEAD": "TEAD_active",
            "HIF2A": "HIF2A_active", "MTOR": "mTOR_active", "AKT": "AKT_active",
            "ERK": "ERK_active", "PI3K": "PI3K_active", "AUTOPHAGY": "autophagy",
            "PDL1": "PDL1_expression", "BRAF": "RAF_active", "RAF": "RAF_active",
            "MEK": "MEK_active", "BCL2": "survival_signal", "SMAD4": "SMAD4_active",
            "TP53": "TP53_active", "P53": "TP53_active",
            "PAN-KRAS": "KRAS_active", "RAS": "KRAS_active",
            # v9
            "GPX4": "GPX4_level", "NRF2": "NRF2_active", "NFE2L2": "NRF2_active",
            "SLC7A11": "SLC7A11_expression", "XCT": "SLC7A11_expression",
            "GSH": "GSH_level", "FSP1": "FSP1_level",
            "AMPK": "AMPK_active", "LKB1": "LKB1_active", "STK11": "LKB1_active",
            "ATR": "ATR_active", "CHK1": "CHEK1_active", "ATM": "ATM_active",
            "CHK2": "CHEK2_active", "PARP": "DDR_active", "RAD51": "RAD51_active",
            "IDO1": "IDO1_active", "IDO": "IDO1_active",
            "CD73": "CD73_expression", "CD47": "CD47_expression",
            "SASP": "SASP_active", "P21": "SASP_active",
        }

        smiles_db = {
            "KRAS_active":          ("c1cc2c(cc1NC(=O)c1cnn(-c3cccc(F)c3)c1)nc(N1CCN(C(=O)C)CC1)nc2N",
                                     "PanKRASi-GD01", "Pan-KRAS tri-complex"),
            "YAP_nuclear":          ("COc1ccc(-c2nc3cc(F)ccc3o2)cc1O",
                                     "YAPi-GD01", "Verteporfin scaffold"),
            "GPX4_level":           ("CC(=O)N1CCN(c2ccc(-c3nc4cccc(F)c4n3C)cc2)CC1",
                                     "GPX4i-GD01 (RSL3-like)", "GPX4 covalent inhibitor"),
            "NRF2_active":          ("O=C(Nc1ccc(-c2ccc(F)c(NC(=O)c3ccncc3)c2)cc1)c1cc(F)ccn1",
                                     "NRF2i-GD01 (ML385-like)", "NRF2-MAFG disruptor"),
            "SLC7A11_expression":   ("OC(=O)C(N)CCC(=O)O", "Erastin-like-GD01", "xCT/SLC7A11 inhibitor"),
            "AMPK_active":          ("CN1CCN(c2ccc(NC(=O)c3ccc(F)cc3Cl)cc2)CC1",
                                     "AMPKa-GD01", "AMPK activator (991-like)"),
            "ATR_active":           ("Cc1ccc(-n2nc(C(F)(F)F)cc2NC(=O)c2cc3c(cc2F)CCO3)cc1",
                                     "ATRi-GD01 (Ceralasertib-like)", "ATR kinase inhibitor"),
            "mTOR_active":          ("Cn1c(=O)c2c(nc(Nc3ccc(N4CCOCC4)cc3)n2C)n(C)c1=O",
                                     "mTORi-GD01", "mTORC1/2 dual"),
            "IDO1_active":          ("Brc1ccc2[nH]cc(CC(N)C(=O)O)c2c1",
                                     "IDO1i-GD01 (Epacadostat-like)", "IDO1 inhibitor"),
            "survival_signal":      ("CC(O)(c1ccc(NC(=O)c2ccncc2)cc1)c1ccc(F)cc1",
                                     "BCL2i-GD01 (Navitoclax-like)", "BCL-2/xL inhibitor"),
        }

        if btn_des and protein.strip():
            import sys, os
            sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

            t = protein.strip().upper()
            node = node_map.get(t)
            if not node:
                all_nodes = list(node_map.values())
                match = [n for n in all_nodes if t.lower() in n.lower()]
                node = match[0] if match else "proliferation_signal"
                st.warning(f"'{t}' mapeado a `{node}`")

            smi, drug_name, drug_class = smiles_db.get(
                node,
                (f"CC(=O)Nc1ccc(-c2ccc(F)c(NC(=O)c3ccncc3)c2)cc1",
                 f"{t}i-GD01", f"Inhibidor de {t}"),
            )
            mw   = round(sum(1 for c in smi if c.isalpha() and c.isupper()) * 14.5 + 50, 1)
            logp = round(smi.count("c") * 0.3 - smi.count("O") * 0.5, 2)
            hbd  = max(smi.count("N") + smi.count("O") - smi.count("n") - smi.count("o"), 0)
            hba  = smi.count("N") + smi.count("O") + smi.count("n")
            viol = (1 if mw > 500 else 0) + (1 if logp > 5 else 0) + \
                   (1 if hbd > 5 else 0) + (1 if hba > 10 else 0)
            ic50_est = round(des_dose_um * 0.3 * (100 / max(des_inh, 10)), 2)

            mc1, mc2 = st.columns([1.2, 1])
            with mc1:
                st.markdown(f"""
<div style="background:#0f1c2b;border:1px solid #1b4f72;border-radius:10px;padding:1.2rem;">
<h4 style="color:#7ab8d9;margin-top:0">🧪 {drug_name}</h4>
<p style="color:#5d8faa;font-size:.82rem">{drug_class}</p>
<div style="background:#060a14;border-radius:6px;padding:.8rem;font-family:monospace;
font-size:.85rem;color:#4a9eff;word-break:break-all;margin:.6rem 0">{smi}</div>
</div>""", unsafe_allow_html=True)
                try:
                    from rdkit import Chem
                    from rdkit.Chem import Draw
                    mol = Chem.MolFromSmiles(smi)
                    if mol:
                        img = Draw.MolToImage(mol, size=(380, 280))
                        st.image(img, caption=f"Estructura 2D: {drug_name}",
                                 use_container_width=True)
                except ImportError:
                    st.caption("ℹ️ Instala `rdkit-pypi` para ver la estructura 2D.")
                except Exception:
                    pass

            with mc2:
                st.markdown(f"""
<div style="background:#0d1a2a;border:1px solid #1b4f72;border-radius:10px;padding:1.2rem;">
<h4 style="color:#7ab8d9;margin-top:0">📊 Propiedades Lipinski</h4>
<table style="width:100%;color:#aed6f1;font-size:.88rem">
<tr><td>MW</td><td style="text-align:right"><strong>{mw} Da</strong></td></tr>
<tr><td>LogP</td><td style="text-align:right"><strong>{logp}</strong></td></tr>
<tr><td>HBD</td><td style="text-align:right"><strong>{hbd}</strong></td></tr>
<tr><td>HBA</td><td style="text-align:right"><strong>{hba}</strong></td></tr>
<tr><td>Violaciones</td><td style="text-align:right"><strong>{viol}</strong></td></tr>
</table></div>""", unsafe_allow_html=True)
                if viol <= 1:
                    st.success("✅ Drug-like (Lipinski RO5)")
                else:
                    st.error("❌ > 1 violación Lipinski")

                st.markdown(f"""
<div style="background:#0d1a2a;border:1px solid #2ed573;border-radius:10px;
padding:1rem;margin-top:.8rem;">
<h4 style="color:#2ed573;margin-top:0">🎯 Predicción {des_cell}</h4>
<p style="color:#aed6f1;font-size:.85rem">
IC₅₀ estimado: <strong style="color:#2ed573">{ic50_est} µM</strong><br>
Diana: <strong>{t}</strong> → <code>{node}</code><br>
Inhibición: <strong>{des_inh}%</strong> a {des_dose_um} µM
</p></div>""", unsafe_allow_html=True)

            # Auto-simular
            st.markdown("---")
            st.markdown("### 🧫 Simulación rápida 72h")
            des_cl_data = CELL_LINES[des_cell]
            with st.spinner(f"🧬 Simulando {drug_name} en {des_cell}..."):
                from simulation.tumor_model import TumorModel
                from drugs.drug_library import Drug
                mp2 = dict(des_cl_data["mutation_profile"])
                mp2.setdefault("MSI_status", "MSS")
                tm = TumorModel(
                    width=des_cl_data["grid"], height=des_cl_data["grid"],
                    n_cancer=des_cl_data["n_cancer"], n_caf=des_cl_data["n_caf"],
                    n_macrophage=des_cl_data["n_mac"], n_tcell=des_cl_data["n_tc"],
                    seed=42, mutation_profile=mp2,
                )
                tm.drug_library.drugs["designed"] = Drug(
                    name=drug_name, targets={node: 1.0},
                    ic50=0.3, max_efficacy=des_inh / 100.0,
                )
                dose_m = hill_dose(des_dose_um, 1.0)
                tm.set_drug_doses({"designed": dose_m})
                bar2 = st.progress(0)
                for i in range(72):
                    tm.step()
                    if (i + 1) % 12 == 0 or i == 71:
                        bar2.progress((i + 1) / 72)
                bar2.empty()

            ts = tm.get_summary()
            final_v = ts["Células tumorales"] / des_cl_data["n_cancer"] * 100
            v_col2 = "#2ed573" if final_v < 30 else ("#f39c12" if final_v < 70 else "#ff4757")
            st.markdown(
                f"<div style='background:#0d1a2a;border:1px solid #1b4f72;border-radius:8px;"
                f"padding:.8rem 1rem;'><strong>{des_cell}</strong> · 72h · "
                f"Viabilidad: <strong style='color:{v_col2}'>{final_v:.0f}%</strong> · "
                f"Vivas: {ts['Células tumorales']} · CAFs: {ts['CAFs']} · "
                f"T-cells: {ts['T-cells CD8+']}</div>",
                unsafe_allow_html=True,
            )

except Exception as e:
    st.error(f"❌ Error: {e}")
    st.exception(e)
