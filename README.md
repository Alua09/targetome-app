# 🎯 Targetome | AI Precision Oncology Platform

**Targetome** is an advanced AI-powered precision oncology platform designed to map patient genomic data into 3D cellular networks, locate the primary target gene driving tumor growth, and provide transparent, data-driven clinical insights.

---

## 🚀 Key Features

* **🌐 3D Tumor Cellular Model**: Interactive Plotly network visualization centered around the core hub gene, built dynamically from the patient's actual mutation profile with adjustable network resolution.
* **📈 Multi-Omic & Survival Analytics**: Comprehensive clinical telemetry, Kaplan-Meier survival curve forecasts, mutation type distributions, and copy number variations (CNA).
* **📋 Mutation Telemetry**: Transparent data tables displaying raw and processed patient genetic records directly from TCGA datasets.
* **🤖 AI Clinical Co-Pilot**: Interactive Q&A interface featuring quick clinical query buttons and natural language processing to explain evidence tiers (OncoKB), targeted therapies, and biological mechanisms.

---

## 🔬 Scientific & Analytical Integrity

Targetome is built on a **100% patient-driven architecture**:
* **No Artificial Padding**: Networks and hub identifications rely strictly on unique patient genes and their mutation frequencies found in the uploaded `.tsv` dataset.
* **Dynamic Resolution**: The network resolution slider filters visual graphs based on real patient-specific nodes.
* **Real-Time Telemetry**: Automatic parsing of TCGA datasets ensures accurate, real-time node and mutation counting without cross-session caching issues.

---

## 📦 Installation & Local Run

1. Clone the repository:
   ```bash
   git clone [https://github.com/your-username/targetome-app.git](https://github.com/your-username/targetome-app.git)
   cd targetome-app
