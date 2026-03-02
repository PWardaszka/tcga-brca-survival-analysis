# =========================================================
# TCGA BRCA Survival Analysis - Streamlit App
# Kaplan–Meier Survival Curves
# =========================================================

import os
import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

# ---------------------------------------------------------
# Page configuration
# ---------------------------------------------------------

st.set_page_config(
    page_title="TCGA BRCA Survival Analysis",
    layout="wide"
)

st.title("Breast Cancer Survival Analysis")
st.markdown(
"""
Kaplan–Meier survival analysis based on TCGA clinical data.
Stratify patients by stage, age group or treatment type.
"""
)

# ---------------------------------------------------------
# Load data
# ---------------------------------------------------------

@st.cache_data
def load_data():
    
    BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    DATA_PATH = os.path.join(BASE_DIR, "data", "survival_cleaned.csv")
    
    df = pd.read_csv(DATA_PATH)
    return df

df = load_data()

# ---------------------------------------------------------
# Basic cleaning
# ---------------------------------------------------------

df = df.dropna(subset=["time", "event"])
df["event"] = df["event"].astype(int)

# Create age groups if not already present
if "age_group" not in df.columns:
    df["age_group"] = df["age"].apply(
        lambda x: "<50" if x < 50 else "50-64" if x < 65 else "≥65"
    )

# ---------------------------------------------------------
# Sidebar - variable selection with treatment options
# ---------------------------------------------------------

st.sidebar.header("Stratification")

# Therapy options extracted from treatments to show in sidebar
therapy_options = [
    "Surgery",
    "Chemotherapy",
    "Radiation",
    "Hormone Therapy",
    "Targeted Molecular Therapy"
]

variable = st.sidebar.selectbox(
    "Select variable for survival comparison:",
    ["stage_group", "age_group", "therapy"]
)

# ---------------------------------------------------------
# Prepare dataframe for plotting based on selection
# ---------------------------------------------------------

if variable == "therapy":
    # Select therapy from sidebar
    therapy_option = st.sidebar.selectbox(
        "Select therapy to analyze:",
        therapy_options
    )
    # Create binary column if patient had this therapy or not
    df["had_therapy"] = np.where(
        df["treatment"].str.contains(therapy_option, case=False, na=False),
        1,
        0
    )
    df_plot = df.dropna(subset=["had_therapy", "time", "event"])
    plot_variable = "had_therapy"
else:
    df_plot = df.dropna(subset=[variable])
    plot_variable = variable

# ---------------------------------------------------------
# Kaplan–Meier Plot
# ---------------------------------------------------------

st.subheader(f"Kaplan–Meier Survival by {variable if variable != 'therapy' else therapy_option}")

kmf = KaplanMeierFitter()

fig, ax = plt.subplots(figsize=(8,6))

for group in sorted(df_plot[plot_variable].unique()):
    mask = df_plot[plot_variable] == group
    label = str(group)
    # For therapy variable, replace 0/1 with No / Yes
    if plot_variable == "had_therapy":
        label = "No " + therapy_option if group == 0 else therapy_option

    kmf.fit(
        durations=df_plot.loc[mask, "time"],
        event_observed=df_plot.loc[mask, "event"],
        label=label
    )
    kmf.plot_survival_function(ax=ax)

ax.set_title(f"Survival Curves Stratified by {variable if variable != 'therapy' else therapy_option}")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Survival Probability")
ax.grid(True)

st.pyplot(fig)

# ---------------------------------------------------------
# Log-rank test
# ---------------------------------------------------------

st.subheader("Log-rank Test")

results = multivariate_logrank_test(
    df_plot["time"],
    df_plot[plot_variable],
    df_plot["event"]
)

st.write(results.summary)

# ---------------------------------------------------------
# Median survival
# ---------------------------------------------------------

st.subheader("Median Survival Time (days)")

median_table = []

for group in sorted(df_plot[plot_variable].unique()):
    
    mask = df_plot[plot_variable] == group
    
    kmf.fit(
        durations=df_plot.loc[mask, "time"],
        event_observed=df_plot.loc[mask, "event"]
    )
    
    median_table.append({
        "Group": ("No " + therapy_option if group == 0 else therapy_option) if plot_variable == "had_therapy" else group,
        "Median Survival (days)": kmf.median_survival_time_
    })

st.dataframe(pd.DataFrame(median_table))

# ---------------------------------------------------------
# Footer
# ---------------------------------------------------------

st.markdown("---")
st.markdown("Data source: TCGA Breast Cancer (TCGA-BRCA)")