"""
🧬 PDAC Gemelo Digital v4 — Entry Point
Ejecuta: streamlit run main.py
"""
import streamlit as st
import sys
import time
import os

# ═══════════════════════════════════════
#  CONFIGURAR PÁGINA — PRIMERA LÍNEA
# ═══════════════════════════════════════
st.set_page_config(page_title="PDAC Gemelo Digital 2026", layout="wide", page_icon="🧬",
                   initial_sidebar_state="expanded")

# ═══════════════════════════════════════
#  TÍTULO — INMEDIATO (anti pantalla blanca)
# ═══════════════════════════════════════
st.title("🧬 PDAC In Silico — Gemelo Digital 2026")
st.markdown("**Simulador de células PDAC + diseño de fármacos contra nuevas dianas**")

# ═══════════════════════════════════════
#  PATH setup
# ═══════════════════════════════════════
ROOT = os.path.dirname(os.path.abspath(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

# ═══════════════════════════════════════
#  CSS
# ═══════════════════════════════════════
st.markdown("""<style>
section[data-testid="stSidebar"]{background:linear-gradient(180deg,#0f0f2e,#1a1a4e)}
section[data-testid="stSidebar"] * {color: #e8e8f0 !important;}
section[data-testid="stSidebar"] label {color: #ffffff !important; font-weight: 500;}
section[data-testid="stSidebar"] .stSlider label {color: #d0d0e0 !important;}
section[data-testid="stSidebar"] small, section[data-testid="stSidebar"] .stCaption {color: #a0a0c0 !important;}
section[data-testid="stSidebar"] h1,section[data-testid="stSidebar"] h2,section[data-testid="stSidebar"] h3 {color: #ffffff !important;}
section[data-testid="stSidebar"] .stMarkdown p {color: #d8d8f0 !important;}
section[data-testid="stSidebar"] summary {color: #ffffff !important; font-weight: 600;}
section[data-testid="stSidebar"] .stCheckbox label span {color: #e0e0f0 !important;}
section[data-testid="stSidebar"] .stSelectbox label {color: #ffffff !important;}
#MainMenu,footer,header{visibility:hidden}
section[data-testid="stSidebar"]{min-width:340px !important;}
</style>""", unsafe_allow_html=True)

# ═══════════════════════════════════════
#  SIDEBAR — Presets + Parámetros + Fármacos
# ═══════════════════════════════════════
st.sidebar.header("🔬 Configuración del Gemelo")

# ═══════════════════════════════════════
#  PASO 1: PERFIL DEL PACIENTE
# ═══════════════════════════════════════
st.sidebar.markdown("### Paso 1: Perfil Clínico (Tumor Board)")

# ── PRESETS DE TUMOR ──
TUMOR_PRESETS = {
    "(personalizado)": None,
    "🟡 Clásico / Progenitor (~50%)": {
        'desc': "Subtipo más frecuente. Bien diferenciado, estroma moderado. "
                "KRAS G12D + TP53 + CDKN2A. Responde mejor a quimioterapia.",
        'n_cancer': 30, 'n_caf': 12, 'n_mac': 8, 'n_tc': 10,
        'm_kras': True, 'kras_var': 0,  # G12D
        'm_tp53': True, 'm_cdk': True, 'm_sm4': False,
        'm_yap': False, 'm_brc': False, 'm_ptn': False, 'm_myc': False,
    },
    "🔴 Basal-like / Escamoso (~25%)": {
        'desc': "Agresivo, pobremente diferenciado. EMT activa, alta desmoplasia. "
                "KRAS G12V + TP53 + CDKN2A + SMAD4 lost + YAP amp. Peor pronóstico.",
        'n_cancer': 45, 'n_caf': 25, 'n_mac': 12, 'n_tc': 5,
        'm_kras': True, 'kras_var': 1,  # G12V
        'm_tp53': True, 'm_cdk': True, 'm_sm4': True,
        'm_yap': True, 'm_brc': False, 'm_ptn': False, 'm_myc': True,
    },
    "🔵 Inmunogénico / BRCA+ (~10%)": {
        'desc': "Mayor infiltración inmune. BRCA mutado → sensible a PARPi y platinos. "
                "Potencial candidato a inmunoterapia.",
        'n_cancer': 25, 'n_caf': 8, 'n_mac': 6, 'n_tc': 20,
        'm_kras': True, 'kras_var': 0,  # G12D
        'm_tp53': True, 'm_cdk': True, 'm_sm4': False,
        'm_yap': False, 'm_brc': True, 'm_ptn': False, 'm_myc': False,
    },
    "🟠 Desmoplásico / Estroma-rico (~15%)": {
        'desc': "Dominado por CAFs y ECM densa. KRAS G12R (dependiente de macropinocitosis). "
                "Pobre penetración de fármacos. Responde peor a quimioterapia.",
        'n_cancer': 20, 'n_caf': 35, 'n_mac': 15, 'n_tc': 4,
        'm_kras': True, 'kras_var': 2,  # G12R
        'm_tp53': True, 'm_cdk': True, 'm_sm4': True,
        'm_yap': False, 'm_brc': False, 'm_ptn': True, 'm_myc': False,
    },
}

with st.sidebar.expander("🏥 Modelos de Tumor Predeterminados", expanded=False):
    preset_name = st.selectbox("Tipo de tumor", list(TUMOR_PRESETS.keys()), key="preset_sel")
    if TUMOR_PRESETS[preset_name] is not None:
        st.caption(TUMOR_PRESETS[preset_name]['desc'])
    btn_preset = st.button("📋 Cargar Preset", use_container_width=True,
        disabled=(preset_name == "(personalizado)"))

# Aplicar preset
if btn_preset and TUMOR_PRESETS[preset_name] is not None:
    p = TUMOR_PRESETS[preset_name]
    for k, v in p.items():
        if k in ('desc', 'kras_var'):
            continue
        if k.startswith(('n_', 'm_')):
            st.session_state[k] = v
    # Set KRAS variant by name
    variant_names = ["G12D (~41%)", "G12V (~30%)", "G12R (~16%)",
                     "G12C (~2%)", "Q61H (~2%)", "WT (wild-type)"]
    st.session_state['kras_var'] = variant_names[p['kras_var']]
    st.rerun()

# ── Mutaciones ──
with st.sidebar.expander("🧬 Perfil Mutacional", expanded=False):
    st.caption("Activa/desactiva mutaciones driver del PDAC")
    st.info("⚠️ Los cambios mutacionales rigen la biología desde el inicio. Para aplicarlos debes iniciar simulación.")

    # KRAS con variante específica
    mut_kras = st.checkbox("KRAS mutado (~92%)", value=True, key="m_kras",
        help="Oncogén driver principal del PDAC. Sin él, el tumor colapsa.")
    kras_variant = st.selectbox("Variante KRAS", [
        "G12D (~41%)", "G12V (~30%)", "G12R (~16%)",
        "G12C (~2%)", "Q61H (~2%)", "WT (wild-type)",
    ], key="kras_var", disabled=not mut_kras,
        help="G12D: Mayor dependencia PI3K. G12R: Pobre macropinocitosis.")
    
    kras_var_code = kras_variant.split(" ")[0] if mut_kras else "WT"

    mc1, mc2 = st.columns(2)
    mut_tp53 = mc1.checkbox("TP53 loss (~75%)", value=True, key="m_tp53", help="Pierde apoptosis.")
    mut_cdkn2a = mc1.checkbox("CDKN2A del (~90%)", value=True, key="m_cdk", help="Pierde freno celular.")
    mut_smad4 = mc1.checkbox("SMAD4 loss (~55%)", value=True, key="m_sm4", help="Promueve metástasis.")
    mut_yap = mc2.checkbox("YAP amplif (~30%)", value=False, key="m_yap", help="Resistencia a inhibidores KRAS.")
    mut_brca = mc2.checkbox("BRCA1/2 mut (~6%)", value=False, key="m_brc", help="Sensibilidad a Platino/PARP.")
    mut_pten = mc2.checkbox("PTEN loss (~5%)", value=False, key="m_ptn", help="Hiperactiva vía mTOR directo.")
    mut_myc = mc2.checkbox("MYC amplif (~15%)", value=False, key="m_myc", help="Hiperproliferación constante.")

mutation_profile = {'KRAS': mut_kras, 'KRAS_variant': kras_var_code, 'TP53': mut_tp53, 'CDKN2A': mut_cdkn2a,
                    'SMAD4': mut_smad4, 'YAP_amp': mut_yap, 'BRCA': mut_brca, 'PTEN': mut_pten, 'MYC': mut_myc}

# ═══════════════════════════════════════
#  PASO 2: COMPOSICIÓN DEL CULTIVO
# ═══════════════════════════════════════
st.sidebar.markdown("---")
st.sidebar.markdown("### Paso 2: Microambiente Tumoral")

with st.sidebar.expander("📐 Contadores Celulares", expanded=False):
    st.info("**💡 ¿Por qué modificar la composición?**\n\n"
            "- **CAFs (Fibroblastos)**: Segregan colágeno duro. Si los aumentas, crean un \n'escudo' (Desmoplasia) y la quimio no penetrará.\n\n"
            "- **Macrófagos**: Actúan como basureros. En altos números suprimen a los Linfocitos T.\n\n"
            "- **T-cells (CD8+)**: Soldados del sistema inmune. Si hay demasiados Tregs/CAFs, quedarán bloqueados fuera del tumor.\n\n"
            "- Modificar estos números simula el 'terreno' en el que crecerá el clon cancerígeno inicial.")
    
    grid_size = st.slider("Tamaño de Placa (Grid)", 20, 80, 50, 5, help="Placa más grande = tumor con centro isquémico (hipóxico)")
    sc1, sc2 = st.columns(2)
    n_cancer = sc1.number_input("🔴 Cáncer", 5, 100, 30, 5, key="n_cancer")
    n_caf = sc1.number_input("🟢 CAFs", 0, 50, 15, 5, key="n_caf")
    n_mac = sc2.number_input("🟠 Macrófagos", 0, 30, 10, 5, key="n_mac")
    n_tc = sc2.number_input("🔵 T-cells CD8+", 0, 30, 8, 2, key="n_tc")
    sc3, sc4 = st.columns(2)
    n_treg = sc3.number_input("🟣 Tregs", 0, 20, 0, 1, key="n_treg", help="Impiden que las T-cells ataquen")
    n_mdsc = sc3.number_input("⚫ MDSCs", 0, 20, 0, 1, key="n_mdsc")
    n_nk = sc4.number_input("🟤 NK cells", 0, 15, 0, 1, key="n_nk")

# ═══════════════════════════════════════
#  PASO 3: TRATAMIENTO (FARMACOLOGÍA)
# ═══════════════════════════════════════
st.sidebar.markdown("---")
st.sidebar.markdown("### Paso 3: Estrategia de Tratamiento")

from drugs.treatment_protocols import PROTOCOLS, get_protocol_display_names
protocol_display = get_protocol_display_names()
protocol_choice = st.sidebar.selectbox("📋 Régimen clínico automático", [
    "(Manual — sliders)", *protocol_display.keys()
], help="Los protocolos estándar (ej: FOLFIRINOX) aplican dosis exactas y deshabilitan la entrada manual.")
use_protocol = protocol_choice != "(Manual — sliders)"
protocol_key = protocol_display.get(protocol_choice, None)

disabled_manual = use_protocol

if use_protocol and protocol_key:
    p = PROTOCOLS[protocol_key]
    st.sidebar.success(f"✔️ Protocolo activo. Manual bloqueado.")
    st.sidebar.caption(f"📖 {p.get('reference', '')}")

# ── Quimioterapia ──
with st.sidebar.expander("💊 Quimioterapia", expanded=not use_protocol):
    dose_gem = st.slider("Gemcitabina", 0.0, 1.0, 0.0, 0.05, key="d_gem", disabled=disabled_manual)
    dose_5fu = st.slider("5-Fluorouracilo (5-FU)", 0.0, 1.0, 0.0, 0.05, key="d_5fu", disabled=disabled_manual)
    dose_oxali = st.slider("Oxaliplatino", 0.0, 1.0, 0.0, 0.05, key="d_oxa", disabled=disabled_manual)
    dose_irino = st.slider("Irinotecán", 0.0, 1.0, 0.0, 0.05, key="d_iri", disabled=disabled_manual)
    dose_nabpac = st.slider("Nab-Paclitaxel", 0.0, 1.0, 0.0, 0.05, key="d_nab", disabled=disabled_manual)

# ── Pan-KRAS 2026 ──
with st.sidebar.expander("🧬 Inhibidores KRAS", expanded=not use_protocol):
    pan_kras_drug = st.selectbox("Fármaco Dique", [
        "(ninguno)", "Daraxonrasib (RMC-6236)", "RMC-7977",
        "ERAS-0015", "LY-4066434", "MRTX1133 (G12D)", "Sotorasib (G12C)",
    ], key="pk_drug", disabled=disabled_manual)
    pan_kras_dose = st.slider("Dosis", 0.0, 1.0, 0.0, 0.05, key="pk_dose", disabled=disabled_manual)
    combo_chemo = st.checkbox("Combo Gemcitabina", key="pk_combo", disabled=disabled_manual)

# ── Terapias dirigidas ──
with st.sidebar.expander("🎯 Terapia Dirigida", expanded=False):
    dose_afatinib = st.slider("Afatinib (pan-HER)", 0.0, 1.0, 0.0, 0.05, key="d_afa", disabled=disabled_manual)
    dose_vert = st.slider("Verteporfin (YAPi)", 0.0, 1.0, 0.0, 0.05, key="d_ver", disabled=disabled_manual)
    dose_trametinib = st.slider("Trametinib (MEKi)", 0.0, 1.0, 0.0, 0.05, key="d_tra", disabled=disabled_manual)
    dose_everolimus = st.slider("Everolimus (mTORi)", 0.0, 1.0, 0.0, 0.05, key="d_eve", disabled=disabled_manual)
    dose_olaparib = st.slider("Olaparib (PARPi)", 0.0, 1.0, 0.0, 0.05, key="d_ola", disabled=disabled_manual)

# ── Inmunoterapia ──
with st.sidebar.expander("🛡️ Inmunoterapia", expanded=False):
    dose_pd1 = st.slider("Anti-PD-1 (Pembrolizumab)", 0.0, 1.0, 0.0, 0.05, key="d_pd1", disabled=disabled_manual)
    dose_ctla4 = st.slider("Anti-CTLA-4 (Ipilimumab)", 0.0, 1.0, 0.0, 0.05, key="d_ct4", disabled=disabled_manual)

# ── Experimental ──
with st.sidebar.expander("🧪 Experimental", expanded=False):
    dose_protac = st.slider("PROTAC STAT3 (SD-36)", 0.0, 1.0, 0.0, 0.05, key="d_prot", disabled=disabled_manual)
    dose_belz = st.slider("Belzutifan (HIF-2αi)", 0.0, 1.0, 0.0, 0.05, key="d_bel", disabled=disabled_manual)
    dose_hcq = st.slider("Hidroxicloroquina (AutoI)", 0.0, 1.0, 0.0, 0.05, key="d_hcq", disabled=disabled_manual)

# ── Diana personalizada ──
with st.sidebar.expander("🔬 Inhibidor Sintético Personalizado", expanded=False):
    custom_target = st.selectbox("Proteína a inhibir", [
        "(ninguna)", "KRAS_active", "YAP_nuclear", "TEAD_active",
        "HIF2A_active", "mTOR_active", "AKT_active", "ERK_active",
        "PI3K_active", "autophagy", "macropinocytosis", "PDL1_expression",
    ], disabled=disabled_manual)
    custom_inh = st.slider("% Eficacia bloqueo", 0, 100, 0, 5, disabled=disabled_manual)

st.sidebar.markdown("---")
st.sidebar.subheader("⏱️ Duración de simulación")
time_presets = {
    "🔬 Agudo (3 días)": 72,
    "📅 1 semana": 168,
    "🗓️ 1 mes": 720,
    "📊 3 meses (RECIST)": 2160,
    "📈 6 meses (PFS)": 4320,
}
time_choice = st.sidebar.selectbox("Duración", list(time_presets.keys()),
    help="Tiempo de simulación. RECIST requiere ≥3 meses.")
sim_steps = time_presets[time_choice]
st.sidebar.caption(f"= {sim_steps}h = {sim_steps/24:.0f} días = {sim_steps/720:.1f} meses")

# ═══════════════════════════════════════
#  ESTADO
# ═══════════════════════════════════════
if 'model' not in st.session_state:
    st.session_state.model = None
if 'ran' not in st.session_state:
    st.session_state.ran = False

# ═══════════════════════════════════════
#  TABS
# ═══════════════════════════════════════
tab_sim, tab_clinic, tab_ab, tab_multi, tab_vcf, tab_drug = st.tabs([
    "🧫 Simulación", "🏥 Dashboard Clínico", "⚖️ Comparador A/B",
    "📈 Multi-Brazo", "📂 Importar Paciente", "💡 Diseñador"])

try:
    # ══════════════════════════════════
    #  TAB: SIMULACIÓN
    # ══════════════════════════════════
    with tab_sim:
        st.markdown("")
        _, cc, _ = st.columns([1, 2, 1])
        with cc:
            btn = st.button(f"🚀 INICIAR SIMULACIÓN ({sim_steps//24} días)",
                            type="primary", use_container_width=True)

        if btn:
            import numpy as np
            # Reset completo — seed aleatorio para resultados distintos cada vez
            st.session_state.model = None
            st.session_state.ran = False
            run_seed = int(time.time() * 1000) % (2**31)
            with st.spinner("⏳ Cargando modelo..."):
                from simulation.tumor_model import TumorModel
                from drugs.drug_library import Drug
                model = TumorModel(width=grid_size, height=grid_size,
                    n_cancer=n_cancer, n_caf=n_caf,
                    n_macrophage=n_mac, n_tcell=n_tc,
                    n_treg=n_treg, n_mdsc=n_mdsc, n_nk=n_nk,
                    seed=run_seed, mutation_profile=mutation_profile)

            # Fármacos — todos los sliders
            doses = {
                'gemcitabine': dose_gem,
                '5fu': dose_5fu,
                'oxaliplatin': dose_oxali,
                'irinotecan': dose_irino,
                'nab_paclitaxel': dose_nabpac,
                'afatinib': dose_afatinib,
                'verteporfin': dose_vert,
                'trametinib': dose_trametinib,
                'everolimus': dose_everolimus,
                'olaparib': dose_olaparib,
                'anti_pd1': dose_pd1,
                'anti_ctla4': dose_ctla4,
                'protac_stat3': dose_protac,
                'belzutifan': dose_belz,
                'hydroxychloroquine': dose_hcq,
            }

            # Pan-KRAS selectbox
            pk_map = {
                "Daraxonrasib (RMC-6236)": 'daraxonrasib',
                "RMC-7977": 'rmc7977',
                "ERAS-0015": 'eras0015',
                "LY-4066434": 'ly4066434',
                "MRTX1133 (G12D)": 'mrtx1133',
                "Sotorasib (G12C)": 'sotorasib',
            }
            if pan_kras_drug != "(ninguno)" and pan_kras_dose > 0:
                pk_id = pk_map[pan_kras_drug]
                doses[pk_id] = pan_kras_dose

            if combo_chemo and pan_kras_drug != "(ninguno)":
                doses['gemcitabine'] = max(doses.get('gemcitabine', 0), 0.5)

            # Si se seleccionó protocolo clínico → scheduling automático
            if use_protocol and protocol_key:
                model.set_protocol(protocol_key)
                # El protocolo controla dosis automáticamente en cada step()
            else:
                model.set_drug_doses(doses)

            if custom_target != "(ninguna)" and custom_inh > 0:
                model.drug_library.drugs['custom'] = Drug(
                    name=f"Inhibidor de {custom_target}",
                    targets={custom_target: 1.0}, ic50=0.3, max_efficacy=0.9)
                doses['custom'] = custom_inh / 100.0
                model.set_drug_doses(doses)

            bar = st.progress(0, text="🧬 Simulando...")
            update_freq = max(sim_steps // 50, 5)  # Adaptar frecuencia
            for i in range(sim_steps):
                model.step()
                if (i+1) % update_freq == 0 or i == sim_steps-1:
                    s = model.get_summary()
                    day = (i+1) // 24
                    week = day // 7
                    if sim_steps > 168:
                        time_str = f"Día {day} (Sem {week})"
                    else:
                        time_str = f"Hora {i+1}"
                    bar.progress((i+1)/sim_steps,
                        text=f"{time_str} — Cáncer: {s['Células tumorales']}")
            bar.empty()
            st.session_state.model = model
            st.session_state.ran = True
            st.toast(f"✅ {sim_steps}h ({sim_steps//24}d) completados", icon="🧬")
            st.rerun()

        # ── Resultados ──
        model = st.session_state.model
        if model is not None and st.session_state.ran:
            import numpy as np
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
            import pandas as pd

            s = model.get_summary()
            st.markdown("---")
            c1,c2,c3,c4,c5,c6 = st.columns(6)
            c1.metric("🕐 Paso", s['Paso'])
            c2.metric("🔴 Cáncer", s['Células tumorales'])
            c3.metric("🟢 CAFs", s['CAFs'])
            c4.metric("🟠 Macrófagos", s['Macrófagos M2'])
            c5.metric("🔵 T-cells", s['T-cells CD8+'])
            c6.metric("🫁 O₂", s['O₂ promedio'])

            # Grid
            grid = model.get_grid_state()
            h, w = grid.shape
            img = np.zeros((h, w, 3), dtype=np.uint8)
            for v, c in {0:(10,10,30),1:(255,71,87),2:(87,96,111),
                         3:(46,213,115),4:(255,165,2),5:(30,144,255)}.items():
                img[grid==v] = c

            gc, oc = st.columns(2)
            with gc:
                fig = go.Figure(go.Image(z=img))
                fig.update_layout(title="🔬 Grid Tumoral",
                    margin=dict(l=10,r=10,t=40,b=10), height=450,
                    xaxis=dict(visible=False), yaxis=dict(visible=False))
                st.plotly_chart(fig, use_container_width=True, key="grid")
            with oc:
                fig2 = go.Figure(go.Heatmap(z=model.microenv.oxygen.T,
                    colorscale=[[0,'#1a0a2e'],[.2,'#6b1d9e'],[.4,'#e74c3c'],
                                [.6,'#f39c12'],[.8,'#2ecc71'],[1,'#3498db']],
                    zmin=0,zmax=1))
                fig2.update_layout(title="🫁 Oxígeno",
                    margin=dict(l=10,r=10,t=40,b=10), height=450,
                    xaxis=dict(visible=False), yaxis=dict(visible=False))
                st.plotly_chart(fig2, use_container_width=True, key="oxy")

            # Evolución
            if model.history['step']:
                df = pd.DataFrame(model.history)
                st.markdown("---")
                fig3 = make_subplots(rows=2,cols=2,
                    subplot_titles=["Poblaciones","Señalización","Subtipo %","PD-L1"],
                    vertical_spacing=.15,horizontal_spacing=.1)
                fig3.add_trace(go.Scatter(x=df.step,y=df.cancer_alive,name='Cáncer',
                    line=dict(color='#ff4757',width=2),fill='tozeroy',
                    fillcolor='rgba(255,71,87,.1)'),row=1,col=1)
                fig3.add_trace(go.Scatter(x=df.step,y=df.caf_count,name='CAFs',
                    line=dict(color='#2ed573',width=2)),row=1,col=1)
                fig3.add_trace(go.Scatter(x=df.step,y=df.tcell_count,name='T-cells',
                    line=dict(color='#1e90ff',width=2)),row=1,col=1)
                fig3.add_trace(go.Scatter(x=df.step,y=df.macrophage_count,name='Mac',
                    line=dict(color='#ffa502',width=2)),row=1,col=1)
                fig3.add_trace(go.Scatter(x=df.step,y=df.avg_proliferation,name='Prolif',
                    line=dict(color='#ff6b9d'),showlegend=False),row=1,col=2)
                fig3.add_trace(go.Scatter(x=df.step,y=df.avg_apoptosis,name='Apop',
                    line=dict(color='#a29bfe'),showlegend=False),row=1,col=2)
                fig3.add_trace(go.Scatter(x=df.step,y=df.basal_like_pct,name='Basal%',
                    line=dict(color='#e17055'),showlegend=False),row=2,col=1)
                fig3.add_trace(go.Scatter(x=df.step,y=df.avg_pdl1,name='PD-L1',
                    line=dict(color='#fdcb6e'),showlegend=False),row=2,col=2)
                fig3.update_layout(height=500,margin=dict(l=40,r=20,t=40,b=40),
                    legend=dict(orientation="h",y=-.15,x=.5,xanchor="center"))
                st.plotly_chart(fig3,use_container_width=True,key="evo")

                # Leyenda + Export + INFORME
                le, ri = st.columns(2)
                with le:
                    st.markdown("🔴 Cáncer · ⚫ Muerta · 🟢 CAF · 🟠 Mac · 🔵 T-cell")
                    st.markdown("**TCGA:** KRAS ~92% · TP53 ~75% · CDKN2A ~90% · SMAD4 ~55%")
                with ri:
                    st.download_button("📥 CSV", df.to_csv(index=False),
                        "pdac.csv","text/csv",use_container_width=True)

                # ── INFORME CIENTÍFICO Y VÍAS MOLECULARES ──
                st.markdown("---")
                st.header("📄 INFORME CIENTÍFICO Y MOLECULAR")
                
                try:
                    from simulation.report_generator import generate_report, get_pathway_deltas, generate_pathway_network_fig
                    
                    # 1. Infografía de Red Molecular
                    st.markdown("### 1. Estado de la Red Celular (Post-Tratamiento)")
                    st.markdown("Visualización estructural de las vías de señalización de la célula tumoral al finalizar la simulación.")
                    fig_net = generate_pathway_network_fig(model)
                    if fig_net:
                        st.plotly_chart(fig_net, use_container_width=True, key="net_chart")
                        st.info("💡 **Guía de colores:** Nodos en **rojo** indican vías de señalización hiperactivas (mecanismos de escape celular o sobrevida). Nodos en **azul** denotan vías exitosamente bloqueadas o en reposo.")
                    
                    st.markdown("---")
                    
                    # 2. Gráfico de Deltas (Cultivo Celular in silico)
                    st.markdown("### 2. Análisis de Cambio Funcional (Delta T0 vs Tf)")
                    deltas = get_pathway_deltas(model)
                    if deltas:
                        nodes = list(deltas.keys())
                        changes = [v[0] for v in deltas.values()]
                        
                        # Limpiar nombres para UI
                        display_nodes = [n.replace('_active', '').replace('_expression', '').replace('_nuclear', '').replace('_signal', '') for n in nodes]
                        colors = ['#2ecc71' if c < 0 else '#e74c3c' for c in changes]
                        
                        fig_delta = go.Figure()
                        fig_delta.add_trace(go.Bar(
                            x=changes, y=display_nodes,
                            orientation='h',
                            marker_color=colors,
                            text=[f"{c:+.1f}%" for c in changes],
                            textposition='auto'
                        ))
                        fig_delta.update_layout(
                            title="Desviación de actividades respecto al inicio de la simulación",
                            xaxis_title="Delta % (T0 -> Tf)",
                            yaxis_title="Nodo Molecular",
                            height=400, margin=dict(l=20,r=20,t=40,b=20),
                            shapes=[dict(type='line', y0=-1, y1=len(nodes), x0=0, x1=0, line=dict(color='white', width=1))]
                        )
                        st.plotly_chart(fig_delta, use_container_width=True, key="delta_chart")
                        st.caption("🟢 **Barras Verdes (Izquierda):** Vías exitosamente inhibidas relativas al inicio. | 🔴 **Barras Rojas (Derecha):** Vías con hiperactivación adaptativa al tratamiento.")
                        
                    st.markdown("---")
                    
                    # 3. Renderizado de Reporte Extendido
                    st.markdown("### 3. Reporte Diagnóstico Detallado")
                    report = generate_report(model)
                    st.markdown(report)
                    
                    try:
                        import markdown as md
                        from fpdf import FPDF
                        # Limpiar caracteres Unicode no soportados por Helvetica
                        rep_clean = report.replace('⚠', '[⚠ ALERTA]').replace('✓', '[✓ OK]').replace('►', '>').replace('•', '-')
                        rep_clean = rep_clean.replace('↓', '(Baja)').replace('↑', '(Sube)').replace('«', '"').replace('»', '"')
                        rep_clean = rep_clean.replace('—', '-').replace('→', '->').replace('←', '<-')
                        rep_clean = rep_clean.replace('α', 'alpha').replace('β', 'beta').replace('γ', 'gamma').replace('κ', 'kappa')
                        rep_clean = rep_clean.replace('μ', 'u').replace('³', '^3').replace('₅', '5').replace('₀', '0')
                        html = md.markdown(rep_clean)
                        
                        pdf = FPDF()
                        pdf.add_page()
                        pdf.set_font("helvetica", size=10)
                        pdf.write_html(html)
                        pdf_bytes = pdf.output()
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            st.download_button("📄 Descargar PDF Científico",
                                data=pdf_bytes, file_name="pdac_informe_cientifico.pdf",
                                mime="application/pdf", use_container_width=True)
                        with col2:
                            st.download_button("📝 Descargar Markdown",
                                report, "pdac_informe.md", "text/markdown",
                                use_container_width=True)
                    except Exception as re:
                        st.warning(f"Aviso PDF: {re}")
                        st.download_button("📄 Descargar Informe (TXT)",
                            report, "pdac_informe.txt", "text/plain",
                            use_container_width=True)
                except Exception as ex:
                    st.error(f"Error crítico en informe: {ex}")
        else:
            st.markdown("")
            st.info("👈 Configura parámetros en la sidebar y pulsa **🚀 INICIAR SIMULACIÓN**")

    # ══════════════════════════════════
    #  TAB: DASHBOARD CLÍNICO
    # ══════════════════════════════════
    with tab_clinic:
        st.markdown("## 🏥 Dashboard Clínico")

        if st.session_state.ran and st.session_state.model is not None:
            model = st.session_state.model
            h = model.history

            from simulation.clinical_endpoints import (
                classify_recist, get_ca199_trend, estimate_ecog,
                estimate_metastasis_risk, ECOG_DESCRIPTIONS
            )

            # ── RECIST ──
            if h['cancer_alive']:
                baseline = h['cancer_alive'][0]
                current = h['cancer_alive'][-1]
                cat, label, color, emoji = classify_recist(baseline, current)
                change_pct = ((current - baseline) / max(baseline, 1)) * 100

                st.markdown("### 📊 Respuesta Tumoral (RECIST 1.1)")
                rc1, rc2, rc3, rc4 = st.columns(4)
                rc1.metric("Clasificación", f"{emoji} {cat}")
                rc2.metric("Cambio tumoral", f"{change_pct:+.0f}%")
                rc3.metric("Baseline", f"{baseline} cél.")
                rc4.metric("Actual", f"{current} cél.")

                st.markdown(f"**{label}** — {cat}")
                st.markdown("---")

            # ── CA 19-9 ──
            if h.get('ca199') and len(h['ca199']) > 1:
                st.markdown("### 🔬 CA 19-9 (Biomarcador Sérico)")
                import plotly.graph_objects as go

                days = [s/24 for s in h['step']]
                fig_ca = go.Figure()
                fig_ca.add_trace(go.Scatter(
                    x=days, y=h['ca199'], mode='lines',
                    name='CA 19-9', line=dict(color='#ff6d00', width=2)))
                fig_ca.add_hline(y=37, line_dash="dash", line_color="green",
                    annotation_text="Límite normal (37 U/mL)")
                fig_ca.update_layout(
                    xaxis_title="Días", yaxis_title="CA 19-9 (U/mL)",
                    height=300, margin=dict(l=40,r=20,t=20,b=40))
                st.plotly_chart(fig_ca, use_container_width=True, key="ca199_chart")
                st.markdown("---")

            # ── ECOG + Metástasis ──
            if h['cancer_alive']:
                st.markdown("### 🏥 Estado del Paciente")
                ec1, ec2 = st.columns(2)

                ecog = estimate_ecog(current, baseline)
                with ec1:
                    st.metric("ECOG Performance Status", f"{ecog}/4")
                    st.caption(ECOG_DESCRIPTIONS.get(ecog, ""))

                met_risk = estimate_metastasis_risk(h, model.mutation_profile)
                with ec2:
                    risk_color = "🟢" if met_risk < 30 else "🟡" if met_risk < 60 else "🔴"
                    st.metric("Riesgo de Metástasis", f"{risk_color} {met_risk:.0f}%")
                    st.caption("Basado en EMT, invasión, duración, SMAD4")

                st.markdown("---")

            # ── Protocolo activo ──
            if model.protocol_name:
                from drugs.treatment_protocols import PROTOCOLS
                p = PROTOCOLS.get(model.protocol_name, {})
                st.markdown("### 💊 Protocolo Activo")
                st.markdown(f"**{p.get('name', model.protocol_name)}**")
                st.caption(p.get('description', ''))
                st.caption(f"📖 {p.get('reference', '')}")

                # Timeline visual
                if h['step'] and p.get('drugs'):
                    cycle_days = p['cycle_days']
                    total_days = h['step'][-1] / 24
                    st.markdown(f"**Ciclos completados:** {int(total_days / cycle_days)} "
                                f"de {p.get('max_cycles', '∞')}")
                    drug_list = ", ".join(p['drugs'].keys())
                    st.markdown(f"**Fármacos:** {drug_list}")
        else:
            st.info("Ejecuta una simulación primero para ver el dashboard clínico.")

    # ══════════════════════════════════
    #  TAB: COMPARADOR A/B
    # ══════════════════════════════════
    with tab_ab:
        st.markdown("## ⚖️ Comparador A/B — Dos Tratamientos")
        st.markdown("Compara dos regímenes terapéuticos lado a lado con el **mismo tumor**.")
        st.markdown("---")

        ab1, ab2 = st.columns(2)
        with ab1:
            st.markdown("### 🅰️ Brazo A")
            ab_drug_a = st.selectbox("Fármaco A", [
                "gemcitabine", "5fu", "oxaliplatin", "daraxonrasib",
                "mrtx1133", "sotorasib", "afatinib", "verteporfin",
                "olaparib", "trametinib", "anti_pd1", "hydroxychloroquine",
            ], key="ab_a")
            ab_dose_a = st.slider("Dosis A", 0.0, 1.0, 0.5, 0.05, key="abd_a")
        with ab2:
            st.markdown("### 🅱️ Brazo B")
            ab_drug_b = st.selectbox("Fármaco B", [
                "gemcitabine", "5fu", "oxaliplatin", "daraxonrasib",
                "mrtx1133", "sotorasib", "afatinib", "verteporfin",
                "olaparib", "trametinib", "anti_pd1", "hydroxychloroquine",
            ], index=3, key="ab_b")
            ab_dose_b = st.slider("Dosis B", 0.0, 1.0, 0.5, 0.05, key="abd_b")

        ab_steps = st.slider("Pasos A/B", 20, 150, 72, 10, key="ab_steps")
        if st.button("🚀 Comparar A vs B", type="primary", use_container_width=True, key="ab_btn"):
            import time as _time
            with st.spinner("Simulando brazo A..."):
                from simulation.tumor_model import TumorModel
                seed_ab = int(_time.time()*1000)%(2**31)
                mA = TumorModel(width=grid_size, height=grid_size,
                    n_cancer=n_cancer, n_caf=n_caf, n_macrophage=n_mac,
                    n_tcell=n_tc, n_treg=n_treg, n_mdsc=n_mdsc, n_nk=n_nk,
                    seed=seed_ab, mutation_profile=mutation_profile)
                mA.set_drug_doses({ab_drug_a: ab_dose_a})
                for _ in range(ab_steps):
                    mA.step()
            with st.spinner("Simulando brazo B..."):
                mB = TumorModel(width=grid_size, height=grid_size,
                    n_cancer=n_cancer, n_caf=n_caf, n_macrophage=n_mac,
                    n_tcell=n_tc, n_treg=n_treg, n_mdsc=n_mdsc, n_nk=n_nk,
                    seed=seed_ab, mutation_profile=mutation_profile)
                mB.set_drug_doses({ab_drug_b: ab_dose_b})
                for _ in range(ab_steps):
                    mB.step()

            # Resultados lado a lado
            import pandas as pd
            from plotly.subplots import make_subplots
            import plotly.graph_objects as go

            r1, r2 = st.columns(2)
            sA = mA.get_summary()
            sB = mB.get_summary()
            with r1:
                st.metric("🅰️ Cáncer final", sA['Células tumorales'])
            with r2:
                st.metric("🅱️ Cáncer final", sB['Células tumorales'])

            dfA = pd.DataFrame(mA.history)
            dfB = pd.DataFrame(mB.history)
            fig = make_subplots(rows=1, cols=2, subplot_titles=["🅰️ "+ab_drug_a, "🅱️ "+ab_drug_b])
            fig.add_trace(go.Scatter(x=dfA.step, y=dfA.cancer_alive, name='Cáncer A',
                line=dict(color='#ff4757')), row=1, col=1)
            fig.add_trace(go.Scatter(x=dfB.step, y=dfB.cancer_alive, name='Cáncer B',
                line=dict(color='#1e90ff')), row=1, col=2)
            fig.update_layout(height=350, margin=dict(l=40,r=20,t=40,b=40))
            st.plotly_chart(fig, use_container_width=True, key="ab_chart")

            # Tabla resumen
            st.markdown("### Resumen comparativo")
            comp = pd.DataFrame({
                'Métrica': ['Células finales', 'Células resistentes', 'Exhaustion T-cell',
                            'PD-L1 medio', 'Basal-like %'],
                '🅰️ '+ab_drug_a: [
                    sA['Células tumorales'], mA.history['resistant_count'][-1],
                    f"{mA.history['tcell_exhaustion'][-1]:.2f}",
                    f"{mA.history['avg_pdl1'][-1]:.3f}",
                    f"{mA.history['basal_like_pct'][-1]:.0f}%",
                ],
                '🅱️ '+ab_drug_b: [
                    sB['Células tumorales'], mB.history['resistant_count'][-1],
                    f"{mB.history['tcell_exhaustion'][-1]:.2f}",
                    f"{mB.history['avg_pdl1'][-1]:.3f}",
                    f"{mB.history['basal_like_pct'][-1]:.0f}%",
                ],
            })
            st.dataframe(comp, use_container_width=True, hide_index=True)

    # ══════════════════════════════════
    #  TAB: ENSAYO MULTI-BRAZO
    # ══════════════════════════════════
    with tab_multi:
        st.markdown("## 📈 Ensayo Clínico Virtual (Multi-Brazo)")
        st.markdown("Simula múltiples regímenes terapéuticos simultáneamente sobre el mismo tumor para encontrar la mejor opción.")
        st.markdown("---")

        multi_steps = st.slider("Duración del Ensayo (horas)", 24, 240, 96, 24, key="multi_steps")
        
        arms_to_run = st.multiselect("Selecciona los brazos a simular:",
            ["Control (Sin Tratar)", "FOLFIRINOX (Quimio)", "Gemcitabina + Nab-Paclitaxel", 
             "Daraxonrasib (KRASi)", "Daraxonrasib + Afatinib (Dual)", "Daraxonrasib + Afatinib + SD-36 (Triple)",
             "Olaparib (PARPi)", "RMC-4630 (SHP2i) + Daraxonrasib", "Palbociclib (CDK4/6i)", "Immunoterapia (anti-PD1)"],
            default=["Control (Sin Tratar)", "FOLFIRINOX (Quimio)", "Daraxonrasib (KRASi)", "Daraxonrasib + Afatinib (Dual)"]
        )
        
        if st.button("🚀 Iniciar Ensayo Multi-Brazo", type="primary", use_container_width=True, key="multi_btn"):
            if not arms_to_run:
                st.warning("Selecciona al menos un brazo.")
            else:
                import time as _time
                import pandas as pd
                import plotly.graph_objects as go
                from simulation.tumor_model import TumorModel
                
                seed_multi = int(_time.time()*1000)%(2**31)
                results_df = {}
                progress_text = st.empty()
                progress_bar = st.progress(0)
                
                arm_configs = {
                    "Control (Sin Tratar)": {},
                    "FOLFIRINOX (Quimio)": {'5fu': 0.8, 'oxaliplatin': 0.8, 'irinotecan': 0.8},
                    "Gemcitabina + Nab-Paclitaxel": {'gemcitabine': 0.8, 'nab_paclitaxel': 0.8},
                    "Daraxonrasib (KRASi)": {'daraxonrasib': 0.8},
                    "Daraxonrasib + Afatinib (Dual)": {'daraxonrasib': 0.8, 'afatinib': 0.8},
                    "Daraxonrasib + Afatinib + SD-36 (Triple)": {'daraxonrasib': 0.8, 'afatinib': 0.8, 'protac_stat3': 0.8},
                    "Olaparib (PARPi)": {'olaparib': 0.8},
                    "RMC-4630 (SHP2i) + Daraxonrasib": {'rmc4630': 0.8, 'daraxonrasib': 0.8},
                    "Palbociclib (CDK4/6i)": {'palbociclib': 0.8},
                    "Immunoterapia (anti-PD1)": {'anti_pd1': 0.8}
                }
                
                for i, arm in enumerate(arms_to_run):
                    progress_text.text(f"Simulando brazo: {arm} ({i+1}/{len(arms_to_run)})")
                    m = TumorModel(width=grid_size, height=grid_size,
                        n_cancer=n_cancer, n_caf=n_caf, n_macrophage=n_mac,
                        n_tcell=n_tc, n_treg=n_treg, n_mdsc=n_mdsc, n_nk=n_nk,
                        seed=seed_multi, mutation_profile=mutation_profile)
                    
                    m.set_drug_doses(arm_configs.get(arm, {}))
                    for _ in range(multi_steps):
                        m.step()
                    
                    results_df[arm] = pd.DataFrame(m.history)
                    progress_bar.progress((i + 1) / len(arms_to_run))
                
                progress_text.empty()
                progress_bar.empty()
                
                # Plot
                fig = go.Figure()
                colors = ['#2c3e50', '#e74c3c', '#3498db', '#f1c40f', '#2ecc71', '#9b59b6', '#e67e22', '#1abc9c', '#34495e', '#ff7979']
                
                for i, arm in enumerate(arms_to_run):
                    df = results_df[arm]
                    fig.add_trace(go.Scatter(x=df.step, y=df.cancer_alive, mode='lines', name=arm,
                        line=dict(width=3, color=colors[i % len(colors)])))
                
                fig.update_layout(
                    title="Carga Tumoral (Células Cancerosas) a lo largo del tiempo",
                    xaxis_title="Horas", yaxis_title="Nº de Células Tumorales Viables",
                    height=500, margin=dict(l=40,r=20,t=40,b=20),
                    legend=dict(orientation="h", y=-0.2, x=0.5, xanchor="center")
                )
                st.plotly_chart(fig, use_container_width=True, key="multi_chart")
                
                # Tabla
                st.markdown("### Resultados Finales")
                metrics = []
                for arm in arms_to_run:
                    df = results_df[arm]
                    c_ini = df.cancer_alive.iloc[0]
                    c_fin = df.cancer_alive.iloc[-1]
                    pct = (c_fin - c_ini) / max(1, c_ini) * 100
                    metrics.append({
                        "Brazo": arm,
                        "Células Finales": int(c_fin),
                        "Cambio %": f"{pct:+.1f}%",
                        "Resistentes": int(df.resistant_count.iloc[-1]),
                        "PD-L1": f"{df.avg_pdl1.iloc[-1]:.2f}",
                        "Basal %": f"{df.basal_like_pct.iloc[-1]:.0f}%"
                    })
                
                import pandas as pd
                st.dataframe(pd.DataFrame(metrics), use_container_width=True, hide_index=True)

    # ══════════════════════════════════
    #  TAB: IMPORTAR PACIENTE
    # ══════════════════════════════════
    with tab_vcf:
        st.markdown("## 📂 Importar Perfil Molecular del Paciente")
        st.markdown("""
Pega aquí el resultado de un **panel molecular** (FoundationOne, Tempus, MSK-IMPACT,
o similar). El sistema reconocerá los genes y configurará las mutaciones automáticamente.

**Genes reconocidos:** KRAS, TP53, CDKN2A, SMAD4, BRCA1, BRCA2, PTEN, MYC, YAP1
        """)
        st.markdown("---")

        vcf_text = st.text_area("Pega el resultado del panel molecular",
            placeholder="Ejemplo:\nKRAS p.G12V\nTP53 p.R175H (LOF)\nCDKN2A homozygous deletion\nSMAD4 loss\nBRCA2 c.5946delT",
            height=200, key="vcf_input")

        if st.button("🔬 Analizar y Configurar", key="vcf_btn", use_container_width=True):
            if vcf_text.strip():
                text = vcf_text.upper()
                detected = {}

                # KRAS + variante
                if 'KRAS' in text:
                    detected['KRAS'] = True
                    if 'G12D' in text: detected['KRAS_variant'] = 'G12D'
                    elif 'G12V' in text: detected['KRAS_variant'] = 'G12V'
                    elif 'G12R' in text: detected['KRAS_variant'] = 'G12R'
                    elif 'G12C' in text: detected['KRAS_variant'] = 'G12C'
                    elif 'Q61' in text: detected['KRAS_variant'] = 'Q61H'
                    else: detected['KRAS_variant'] = 'G12D'
                else:
                    detected['KRAS'] = False
                    detected['KRAS_variant'] = 'WT'

                detected['TP53'] = 'TP53' in text
                detected['CDKN2A'] = 'CDKN2A' in text or 'P16' in text
                detected['SMAD4'] = 'SMAD4' in text or 'DPC4' in text
                detected['BRCA'] = 'BRCA1' in text or 'BRCA2' in text
                detected['PTEN'] = 'PTEN' in text
                detected['MYC'] = 'MYC' in text and 'MYCN' not in text
                detected['YAP_amp'] = 'YAP' in text

                st.success("✅ Perfil detectado:")
                for gene, status in detected.items():
                    if gene == 'KRAS_variant':
                        continue
                    emoji = "🔴" if status else "⚪"
                    extra = f" ({detected['KRAS_variant']})" if gene == 'KRAS' and status else ""
                    st.write(f"  {emoji} **{gene}**: {'Mutado' if status else 'WT'}{extra}")

                st.info("👆 Ahora ve a la pestaña **🧫 Simulación**, ajusta las mutaciones "
                        "en la sidebar según este perfil, y pulsa Simular.")
            else:
                st.warning("Pega al menos un gen con su alteración.")

    # ══════════════════════════════════
    #  TAB: DISEÑADOR
    # ══════════════════════════════════
    with tab_drug:
        st.markdown("## 💡 Diseñador de Fármacos — Nuevas Dianas")
        st.markdown("---")

        dc1, dc2 = st.columns(2)
        with dc1:
            protein = st.text_input("🧬 Proteína diana",
                placeholder="YAP, TEAD, CREB, MYC, Pan-KRAS...")
            des_inh = st.slider("% Inhibición", 10, 100, 80, 5, key="di")
            btn_des = st.button("🧪 Generar Molécula + Simular",
                type="primary", use_container_width=True)
        with dc2:
            st.markdown("**Dianas conocidas:**")
            st.code("KRAS, YAP, TEAD, HIF2A, mTOR, AKT,\nERK, PI3K, PDL1, Pan-KRAS,\nCREB, MYC, BRAF, MEK, BCL2")

        if btn_des and protein.strip():
            target = protein.strip().upper()
            node_map = {
                'KRAS':'KRAS_active','YAP':'YAP_nuclear','TEAD':'TEAD_active',
                'HIF2A':'HIF2A_active','MTOR':'mTOR_active','AKT':'AKT_active',
                'ERK':'ERK_active','PI3K':'PI3K_active','PDL1':'PDL1_expression',
                'BRAF':'RAF_active','RAF':'RAF_active','MEK':'MEK_active',
                'CREB':'survival_signal','MYC':'proliferation_signal',
                'BCL2':'survival_signal','PAN-KRAS':'KRAS_active',
                'RAS':'KRAS_active','PANKRAS':'KRAS_active',
            }
            node = node_map.get(target, 'proliferation_signal')
            smi_db = {
                'KRAS_active': 'c1cc2c(cc1NC(=O)c1cnn(-c3cccc(F)c3)c1)nc(N1CCN(C(=O)C)CC1)nc2N',
                'YAP_nuclear': 'COc1ccc(-c2nc3cc(F)ccc3o2)cc1O',
                'TEAD_active': 'CC(=O)Nc1ccc(-c2ccc(F)c(NC(=O)c3ccncc3)c2)cc1',
                'HIF2A_active': 'Oc1ccc(-c2nc(-c3ccc(F)cc3)no2)c(O)c1',
                'mTOR_active': 'Cn1c(=O)c2c(nc(Nc3ccc(N4CCOCC4)cc3)n2C)n(C)c1=O',
            }
            smi = smi_db.get(node, 'CC(=O)Nc1ccc(-c2ccc(F)c(NC(=O)c3ccncc3)c2)cc1')
            name = f"{target}i-GD01"

            st.success(f"✅ **{target}** → `{node}` — {des_inh}% inhibición")
            st.code(f"SMILES: {smi}", language=None)

            # RDKit (opcional)
            try:
                from rdkit import Chem
                from rdkit.Chem import Draw
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    img = Draw.MolToImage(mol, size=(400, 300))
                    st.image(img, caption=f"Estructura: {name}")
            except ImportError:
                st.caption("ℹ️ `pip install rdkit-pypi` para ver estructura 2D")
            except Exception:
                pass

            # Simular
            st.markdown("---")
            with st.spinner(f"🧬 Simulando con {name}..."):
                from simulation.tumor_model import TumorModel
                from drugs.drug_library import Drug
                tm = TumorModel(width=grid_size,height=grid_size,
                    n_cancer=n_cancer,n_caf=n_caf,n_macrophage=n_mac,n_tcell=n_tc,
                    n_treg=n_treg,n_mdsc=n_mdsc,n_nk=n_nk,
                    seed=int(time.time()*1000)%(2**31),mutation_profile=mutation_profile)
                tm.drug_library.drugs['designed'] = Drug(
                    name=name,targets={node:1.0},ic50=0.3,max_efficacy=des_inh/100.0)
                tm.set_drug_doses({'designed':des_inh/100.0})
                bar = st.progress(0)
                for i in range(sim_steps):
                    tm.step()
                    if (i+1)%5==0: bar.progress((i+1)/sim_steps)
                bar.empty()
            ts = tm.get_summary()
            st.session_state.model = tm
            st.session_state.ran = True
            st.success(f"Cáncer: **{ts['Células tumorales']}** | CAFs: {ts['CAFs']} | "
                       f"T-cells: {ts['T-cells CD8+']} | O₂: {ts['O₂ promedio']}")
            st.info("👆 Ve a **🧫 Simulación** para gráficos completos")

except Exception as e:
    st.error(f"❌ Error: {e}")
    st.exception(e)
