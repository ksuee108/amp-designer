import streamlit as st
import pandas as pd
import numpy as np
from design import  algorithms_setup, plot_pareto_fronts_many, plot_pareto_fronts_multi, amino_acid_percentage  # 需確認這兩個在 design.py 中實作
import os

user_home = os.path.expanduser("~")

st.set_page_config(page_title="Antimicrobial peptides design by AI", layout="wide")
st.title("🧬 Antimicrobial peptides design by AI")
st.text("This app allows users to explore and design antimicrobial peptides by multi-objective optimization.")
main_tab1, main_tab2, main_tab3 = st.tabs(["Home", "About this App", "How to Design"])

# --- 固定底部 Footer ---
footer = """
<style>
.footer-text {
    position: fixed;
    left: 0;
    bottom: 0;
    width: 100%;
    background: rgb(240,240,240);
    color: black;
    text-align: center;
    padding: 10px 0;
    font-size: 14px;
    border-top: 1px solid #ccc;
    z-index: 1000;
}
</style>

<div class="footer-text">

    🚀 AMP Design App © 2025 
    Citation: Yang C-H, Chen Y-L, Cheung T-H, Chuang L-Y. Multi-Objective Optimization Accelerates the De Novo Design of Antimicrobial Peptide for Staphylococcus aureus. International Journal of Molecular Sciences. 2024; 25(24):13688. https://doi.org/10.3390/ijms252413688
</div>
"""
with main_tab1:
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

    # 可選：清除按鈕
    if st.button("🗑️ Clear All Constraints"):
        st.session_state.constraints = []
        st.info("All constraints have been cleared.")# show all the constraints

    if st.session_state.constraints:
        st.subheader("Current Constraints")
        df_constraints = pd.DataFrame(st.session_state.constraints)
        st.dataframe(df_constraints, use_container_width=True)

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
                df = pd.read_csv(f"dataset\\biopython-{Bacteria[0]}.csv")
                for i in range(1, len(Bacteria)):
                    temp_df = pd.read_csv(f"dataset\\biopython-{Bacteria[i]}.csv")
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
                    path=user_home,
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
            algorithms
        )
        
        if selected_algo:

            if selected_algo:
                results = st.session_state.optimization_results[selected_algo]
                res_dict_flipped = results["res_dict"].copy()
                pareto_df_flipped = results["pareto_df"].copy()
                merged_df_flipped = results["merged_df"].copy()
            
            if "Gravy" in optimization_directions:
                if optimization_directions["Gravy"] == 'hydrophobicity':
                    res_dict_flipped['Gravy'] = -res_dict_flipped['Gravy']
                    pareto_df_flipped['Gravy'] = -pareto_df_flipped['Gravy']
                    merged_df_flipped['Gravy'] = -merged_df_flipped['Gravy']

            results = st.session_state.optimization_results[selected_algo]
            
            tab1, tab2, tab3 = st.tabs(["Objectives", "All Results", "Merged Data"])
            
            with tab1:
                st.write(f"**Objective values for {selected_algo}:**")
                st.dataframe(res_dict_flipped)
            
            with tab2:
                st.write(f"**All optimized results for {selected_algo}:**")
                st.dataframe(pareto_df_flipped)
        
            
            with tab3:
                st.write(f"**Merged data for {selected_algo}:**")
                st.dataframe(merged_df_flipped)

    else:
        st.info("No optimization results cached yet. Run optimization first.")

    for algo in algorithms:
        fasta_path = os.path.join(user_home, f"{algo}.fasta")
        with open(fasta_path, "w") as fasta_file:
            for _, row in merged_df_flipped.iterrows():
                fasta_file.write(f">{row['sequence']}\n{row['sequence']}\n")

        pareto_df_flipped.to_csv(os.path.join(user_home, f"{algo} all optimize result.csv"), index=False)
        st.info(f"Optimize pareto front result saved at file: {user_home}\\{algo} all optimize result.csv")

        merged_df_flipped.to_csv(os.path.join(user_home, f"{algo} optimize result.csv"), index=False)
        st.info(f"Optimize pareto front result saved at file: {user_home}\\{algo} optimize result.csv")

        st.info(f"FASTA saved at file: {fasta_path}")
    st.markdown("---")

    if st.button("📊 Plot Results"):
        with st.spinner("Running optimization... This may take a few minutes."):
            amino_acid_percentage(user_home, algorithms)
            if len(optimization_directions) > 3:
                plot_pareto_fronts_many(user_home, algorithms, optimization_directions)
            else:
                plot_pareto_fronts_multi(user_home, algorithms, optimization_directions)

with main_tab2:
    st.header("About this App")
    st.markdown("""
    ### Overview:
    This Streamlit app is designed to facilitate the de novo design of antimicrobial peptides (AMPs) using multi-objective optimization techniques. Users can select target bacteria, define optimization parameters, and choose from various algorithms to generate peptide sequences with desired physicochemical properties.
    ### Features:
    - **Multi-Bacteria Targeting**: Select from multiple bacteria to tailor peptide designs.
    - **Customizable Optimization**: Choose from various algorithms and define specific objectives and constraints.
    - **Comprehensive Results**: View and download optimized peptide sequences along with their properties.
    - **Visualization Tools**: Generate plots to visualize optimization results and amino acid distributions.
    ### Data Sets:
    The app utilizes precomputed datasets containing physicochemical properties of peptides against various bacteria, including *E. coli*, *S. aureus*, *P. aeruginosa*, and *A. baumannii* from the DBAASP.
    ### Intended Users:
    - Researchers in microbiology and bioinformatics.
    - Pharmaceutical developers focusing on antimicrobial agents.
    - Educators and students in related fields.
    ### Citation:
    If you use this app for your research, please cite the following paper:
    Yang C-H, Chen Y-L, Cheung T-H, Chuang L-Y. Multi-Objective Optimization Accelerates the De Novo Design of Antimicrobial Peptide for Staphylococcus aureus. International Journal of Molecular Sciences. 2024; 25(24):13688. https://doi.org/10.3390/ijms252413688
    """)

with main_tab3:
    st.header("How to Design Antimicrobial Peptides using this App")
    st.markdown("""
    1. **Select Bacteria**: Use the sidebar to choose the bacteria you want to target. *You must select at least one.*
    2. **Set Parameters**: Define population size, peptide length, and number of generations.
    3. **Choose Algorithms**: Select one or more optimization algorithms from the sidebar.
    4. **Select Objectives**: Pick at least two physicochemical properties to optimize.
    5. **Define Constraints**: Optionally, set constraints on specific properties.
    6. **Run Optimization**: Click the "Run Optimization" button to start the process.
    7. **View Results**: After completion, view and download the optimized peptide sequences and their properties.
    8. **Plot Results**: Generate visualizations of the optimization results and amino acid distributions.
    """)
    
st.markdown(footer, unsafe_allow_html=True)