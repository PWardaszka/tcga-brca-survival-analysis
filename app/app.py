# app/app.py

import os
import streamlit as st
import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt

# --- Load cleaned survival data ---
# Determine the path relative to this script
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(BASE_DIR, "..", "data", "clinical.tsv")

df = pd.read_csv(DATA_PATH, sep="\t", low_memory=False)

# --- Prepare survival variables ---
time_death = pd.to_numeric(df.get('demographic.days_to_death'), errors='coerce')
time_followup = pd.to_numeric(df.get('diagnoses.days_to_last_follow_up'), errors='coerce')
time = time_death.fillna(time_followup)

if 'demographic.vital_status' in df.columns:
    event = df['demographic.vital_status'].str.lower().map(lambda x: 1 if x=='dead' else 0)
else:
    event = pd.Series(np.nan, index=df.index)

age = pd.to_numeric(df.get('demographic.age_at_index'), errors='coerce')
stage = df.get('diagnoses.ajcc_pathologic_stage', np.nan)
treatment = df.get('treatments.treatment_type', np.nan)

# Combine into a new DataFrame
df_survival = pd.DataFrame({
    'time': time,
    'event': event,
    'age': age,
    'stage': stage,
    'treatment': treatment
}).dropna(subset=['time','event'])

# --- Streamlit UI ---
st.title("TCGA BRCA Survival Analysis")

# Sidebar filters
st.sidebar.header("Filter options")

filter_type = st.sidebar.radio(
    "Group survival by:",
    ('Stage', 'Treatment')
)

if filter_type == 'Stage':
    group_col = 'stage'
    options = df_survival['stage'].dropna().unique()
else:
    group_col = 'treatment'
    options = df_survival['treatment'].dropna().unique()

selected_group = st.sidebar.multiselect(
    f"Select {group_col}(s) to display",
    options,
    default=options[:3]  # show first 3 by default
)

# --- Plot Kaplan-Meier ---
st.header("Kaplan-Meier Survival Curve")

kmf = KaplanMeierFitter()

plt.figure(figsize=(8,5))

for grp in selected_group:
    mask = df_survival[group_col] == grp
    kmf.fit(df_survival.loc[mask, 'time'], df_survival.loc[mask, 'event'], label=grp)
    kmf.plot_survival_function(ci_show=False)

plt.xlabel("Time (days)")
plt.ylabel("Survival probability")
plt.title("Kaplan-Meier Survival Curve")
plt.grid(True)

st.pyplot(plt)