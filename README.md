# 🧬 PDAC InSilico - Gemelo Digital

Simulador de gemelo digital para **adenocarcinoma ductal pancreático (PDAC)** usando modelado basado en agentes y redes de señalización molecular.

## 🚀 Cómo ejecutar

```bash
# 1. Crear entorno virtual
python -m venv venv
source venv/bin/activate  # En Windows: venv\Scripts\activate

# 2. Instalar dependencias
pip install -r requirements.txt

# 3. Lanzar la aplicación
streamlit run main.py
```

## 📁 Estructura del proyecto

```
pdac_gemelo/
├── main.py                  # Lanzador Streamlit
├── agents/
│   ├── cancer_cell.py       # Célula cancerosa con estado molecular
│   ├── caf.py               # Fibroblastos asociados al cáncer
│   └── immune_cells.py      # Macrófagos M2 y linfocitos T CD8+
├── signaling/
│   └── pdac_network.py      # Red de señalización (KRAS, YAP, PI3K, etc.)
├── microenvironment/
│   └── diffusion.py         # Oxígeno, ECM, citoquinas
├── drugs/
│   └── drug_library.py      # Fármacos reales + mecanismos
├── simulation/
│   └── tumor_model.py       # Modelo principal de Mesa
├── gui/
│   └── app.py               # Interfaz Streamlit completa
├── requirements.txt
└── README.md
```

## 🧪 Características

- **Modelado basado en agentes** (Mesa 3.0+): células cancerosas, CAFs, macrófagos M2, linfocitos T CD8+
- **Señalización molecular**: KRAS G12D, YAP/TAZ-TEAD, PI3K/AKT/mTOR, HIF-2α, TGF-β/SMAD4
- **Microambiente tumoral**: hipoxia, ECM densa, citoquinas
- **Fármacos**: Gemcitabina, MRTX1133 (KRASi), Verteporfin (YAPi), anti-PD-1
- **Subtipos**: classical vs basal-like con plasticity
- **Mutaciones reales**: frecuencias TCGA/CPTAC (KRAS ~92%, TP53 ~75%, SMAD4 ~55%, CDKN2A ~90%)

## 📊 Interfaz

Dashboard interactivo con:
- Sliders para ajustar dosis de fármacos
- Visualización del grid tumoral en 2D
- Gráficos de evolución temporal
- Mapa de hipoxia
- Estadísticas en tiempo real
