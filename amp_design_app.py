import streamlit as st
import pandas as pd
import folium
from streamlit_folium import st_folium
import numpy as np
import random
import time
from design import  algorithms_setup, plot_pareto_fronts_many, plot_pareto_fronts_multi, amino_acid_percentage  # 需確認這兩個在 design.py 中實作

st.set_page_config(page_title="Antimicrobial peptides design by AI", layout="wide")
st.title("🧬 Antimicrobial peptides design by AI")
st.text("This app allows users to explore and design antimicrobial peptides by multi-objective optimization.")

# ----------------------------
# Sidebar inputs
# ----------------------------
Bacteria = st.sidebar.multiselect(
    "Select Bacteria",
    options=["E. coli", "S. aureus", "P. aeruginosa", "A. baumannii"]
)
pop_size = st.sidebar.number_input("Population size", min_value=10, max_value=400, value=80, step=10)
length = st.sidebar.number_input("Peptide sequence length", min_value=10, max_value=20, value=10, step=1)
generations = st.sidebar.number_input("Number of generations", min_value=10, max_value=200, value=100, step=1)

st.sidebar.header("Select Algorithms")
algorithms = st.sidebar.multiselect(
    "Choose optimization algorithms",
    options=[
        'NSGA-II', 'NSGA-III', 'R-NSGA-II', 'R-NSGA-III',
        'U-NSGA-III', 'AGE-MOEA', 'AGE-MOEA-II'
    ]
)
# ----------------------------
# Objective selection
# ----------------------------
st.header("Objectives to optimize")
opt = st.multiselect(
    "Select properties to optimize",
    options=[
        'Gravy', 'Instability Index', 'Aliphatic Index', 'Isoelectric point',
        'Net charge', 'Molecular Weight', 'Charge at pH', 'Aromaticity',
        'Secondary structure fraction Helix', 'Secondary structure fraction Turn',
        'Secondary structure fraction Sheet', 'Boman Index'
    ]
)

if len(opt) < 2:
    st.warning("Please select at least two objectives to optimize.")
    st.stop()


optimization_directions = {}
for i in opt:
    if i != "Gravy" :
        optimization_directions[i] = st.selectbox(
            f"Select minima or maxima for {i}",
            options=['Minimize', 'Maximize'],
            key=i
        )
    else:
        optimization_directions[i] = st.selectbox(
            f"Select optimiz hydrophobicity or hydrophilicity for {i}",
            options=['hydrophobicity', 'hydrophilicity'],
            key=i
        )
st.success("Optimization Directions Set:")
st.json(optimization_directions)
# ----------------------------
# Constraints
# ----------------------------
st.markdown("---")
st.header("Constraints")

options=[
        'Gravy', 'Instability Index', 'Aliphatic Index', 'Isoelectric point',
        'Net charge', 'Molecular Weight', 'Charge at pH', 'Aromaticity',
        'Secondary structure fraction Helix', 'Secondary structure fraction Turn',
        'Secondary structure fraction Sheet', 'Boman Index'
    ]

# 初始化 Session State（避免重新整理後消失）
if "constraints" not in st.session_state:
    st.session_state.constraints = []

# 使用者輸入新約束項
constraint_feature = st.selectbox("Select a physicochemical property:", options)
constraint_type = st.radio(
    "Constraint type:",
    [f"(Maximum limit) ≤ {constraint_feature}", f"{constraint_feature} ≥ (Minimum limit)"]
)
constraint_value = st.number_input(
    f"Enter limit value for {constraint_feature}:",
    value=0.0,
    step=0.1,
    format="%.2f"
)

# 加入按鈕
if st.button("➕ Add Constraint"):
    new_constraint = {
        "Feature": constraint_feature,
        "Type": "max" if "≤" in constraint_type else "min",
        "Value": constraint_value,
    }

    # 避免重複加入相同 feature
    existing = [c for c in st.session_state.constraints if c["Feature"] == constraint_feature]
    if existing:
        st.warning(f"⚠️ {constraint_feature} constraint already exists — updated value.")
        st.session_state.constraints = [
            new_constraint if c["Feature"] == constraint_feature else c
            for c in st.session_state.constraints
        ]
    else:
        st.session_state.constraints.append(new_constraint)
        st.success(f"✅ Added: {constraint_feature} ({new_constraint['Type']} = {constraint_value})")

# 顯示所有約束項
if st.session_state.constraints:
    st.subheader("Current Constraints")
    df_constraints = pd.DataFrame(st.session_state.constraints)
    st.dataframe(df_constraints, use_container_width=True)

    # 可選：清除按鈕
    if st.button("🗑️ Clear All Constraints"):
        st.session_state.constraints = []
        st.info("All constraints have been cleared.")

# 最終輸出（可傳入演算法）
constraint_dict_list = st.session_state.constraints
st.write(constraint_dict_list)

# ----------------------------
# Load data based on selected bacteria
# ----------------------------
st.markdown("---")
if len(Bacteria) < 1:
    st.warning("Please select at least one bacteria from the sidebar.")
    st.stop()
else:
    # Initialize session state for dataframe
    if "loaded_df" not in st.session_state:
        st.session_state.loaded_df = None
    
    # Load data only if not cached
    if st.session_state.loaded_df is None:
        try:
            df = pd.read_csv(f"biopython-{Bacteria[0]}.csv")
            for i in range(1, len(Bacteria)):
                temp_df = pd.read_csv(f"biopython-{Bacteria[i]}.csv")
                df = pd.merge(df, temp_df, on='sequence', suffixes=('', f'_{Bacteria[i]}')).drop_duplicates(subset=['sequence'])
            st.session_state.loaded_df = df
            st.success(f"Loaded data for {', '.join(Bacteria)} with {len(df)} sequences.")
        except FileNotFoundError as e:
            st.error(f"Missing file: {e}")
            st.stop()
    else:
        df = st.session_state.loaded_df
        st.success(f"Using cached data for {', '.join(Bacteria)} with {len(df)} sequences.")

st.subheader("Loaded peptide data preview")
st.dataframe(df.head())

# ----------------------------
# Run optimization
# ----------------------------
st.markdown("---")
if st.button("🚀 Run Optimization"):
    with st.spinner("Running optimization... This may take a few minutes."):
        try:
            # algorithm setup
            setup = algorithms_setup(
                df=df,
                algorithms_list=algorithms,
                pop_size=pop_size,
                generations=generations,
                optimization_directions=optimization_directions,
                length=length,
                opt=opt,
                constraint_dict_list=constraint_dict_list
            )
        
            setup.run_optimization()
            setup.run()
        
            st.success("Optimization completed successfully ✅")

        except Exception as e:
            st.error(f"Error during optimization: {e}")
# ----------------------------
# Display cached results
# ----------------------------
st.subheader("📋 Optimization Results (Cached)")

if "optimization_results" in st.session_state and st.session_state.optimization_results:
    selected_algo = st.selectbox(
        "Select algorithm to view results:",
        list(st.session_state.optimization_results.keys())
    )
    
    if selected_algo:
        results = st.session_state.optimization_results[selected_algo]
        
        tab1, tab2, tab3 = st.tabs(["Objectives", "All Results", "Merged Data"])
        
        with tab1:
            st.write(f"**Objective values for {selected_algo}:**")
            st.dataframe(results['res_dict'])
        
        with tab2:
            st.write(f"**All optimized results for {selected_algo}:**")
            st.dataframe(results['pareto_df'])
        
        with tab3:
            st.write(f"**Merged data for {selected_algo}:**")
            st.dataframe(results['merged_df'])
else:
    st.info("No optimization results cached yet. Run optimization first.")
    
st.markdown("---")
if st.button("📊 Plot Results"):
    amino_acid_percentage(algorithms)
    if len(optimization_directions) > 3:
        plot_pareto_fronts_many(algorithms, optimization_directions)
    else:
        plot_pareto_fronts_multi(algorithms, optimization_directions)