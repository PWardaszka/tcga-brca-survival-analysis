# TCGA-BRCA Survival Analysis App

# Project Goal
The goal of this project is to analyze survival of breast cancer patients (TCGA-BRCA) using Kaplan–Meier survival curves. The app allows stratification of patients by age, disease stage, and treatment type.

# Data Source
Data are obtained from [The Cancer Genome Atlas (TCGA) – Breast Invasive Carcinoma (BRCA)](https://portal.gdc.cancer.gov/).  
The dataset includes clinical information such as survival time, patient status, age, tumor stage, and treatments administered.

# Analysis Description
- Data were cleaned and missing values were handled.
- Age groups (`<50`, `50–64`, `≥65`) and stage groups (`I`, `II`, `III`, `IV`) were created.
- Kaplan–Meier analysis is performed with stratification by:
  - Stage (I, II, III, IV)
  - Age group
  - Treatment type (Surgery, Chemotherapy, Radiation, Hormone Therapy, Targeted Molecular Therapy)
- A log-rank test is conducted to compare survival curves across groups, and median survival times are calculated.

# Results
- Patients with more advanced cancer stages have worse survival outcomes.  
- Age groups show minor differences in survival between patients under 50 and those aged 50–64. Survival decreases noticeably for patients aged 65 and above.  
- Among treatment types, "Targeted Molecular Therapy" is associated with the best survival outcomes.



# How to run
pip install -r requirements.txt

streamlit run app/app_v2.py

