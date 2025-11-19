import streamlit as st
import pandas as pd
import numpy as np
from design import  algorithms_setup, plot_pareto_fronts_many, plot_pareto_fronts_multi, amino_acid_percentage
import os
from BioAnalysis import Bio_analysis
from Bio import SeqIO
from io import StringIO

user_home = os.path.expanduser("~")

st.set_page_config(page_title="Antimicrobial peptides design by AI", layout="wide")
st.title("🧬 Antimicrobial peptides design by AI")
st.text("This app allows users to explore and design antimicrobial peptides by multi-objective optimization.")
main_tab1, main_tab2, main_tab3, main_tab4 = st.tabs(["Home", "About this App", "How to Design", "Related Databases and Prediction Websites"])

footer = """
<style>

main > div {
    padding-bottom: 0px !important;
    padding-top: 0px !important;
}

body {
    margin: 0;
    padding-bottom: 60px; 
}

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
    🚀 AMP Design App © 2025<br>
</div>
"""
st.markdown(footer, unsafe_allow_html=True)

with main_tab1:

    # ----------------------------
    # Sidebar inputs
    # ----------------------------
    all_features = []
    with st.sidebar:
        st.header("Select Bacteria")
        Bacteria = st.multiselect(
            "Choose bacteria",
            options=["E. coli", "S. aureus", "P. aeruginosa", "A. baumannii"]
        )

        st.header("Upload Peptide Sequence")
        uploaded_file = st.file_uploader("Upload FASTA or TXT", type=["txt", "fasta", "fa"])

        seq = None
        bio_analysis = None

        if uploaded_file is not None:
            file_type = uploaded_file.name.lower()

            # -------- FASTA --------
            if file_type.endswith(".fasta") or file_type.endswith(".fa"):
                uploaded_file.seek(0)
                fasta_str = uploaded_file.read().decode("utf-8")
                fasta_io = StringIO(fasta_str)

                records = list(SeqIO.parse(fasta_io, "fasta"))

                if len(records) == 0:
                    st.error("FASTA file is empty or invalid.")
                else:
                    for rec in records:
                        seq = str(rec.seq)

                        try:
                            bio_analysis = Bio_analysis(seq)

                            Gravy = bio_analysis.get_gravy()
                            instability_index = bio_analysis.get_instability_index()
                            Aliphatic_Index = bio_analysis.get_aliphatic_index()
                            Boman_index = bio_analysis.get_boman_index()
                            isoelectric_point = bio_analysis.get_isoelectric_point()
                            net_charge = bio_analysis.get_net_charge()
                            molecular_weight = bio_analysis.get_molecular_weight()
                            charge_at_pH = bio_analysis.get_charge_at_pH()
                            aromaticity = bio_analysis.get_aromaticity()
                            sec_H, sec_T, sec_S = bio_analysis.get_secondary_structure_fraction()
                            amphipathicity = bio_analysis.get_amphipathicity()
                            correlation = bio_analysis.get_auto_correlation()
                            covariance = bio_analysis.get_auto_covariance()
                            hydrophobic_moenet = bio_analysis.get_hydrophobic_moenet()
                            mass = bio_analysis.get_mass()
                            mz = bio_analysis.get_mz()
                            SequenceLength = bio_analysis.get_sequenceLength()

                            features = {
                                "Sequence": seq,
                                'Length': SequenceLength,
                                'Gravy': Gravy,
                                'Instability Index': instability_index,
                                'Aliphatic Index': Aliphatic_Index,
                                'Isoelectric point': isoelectric_point,
                                'Net charge': net_charge, 
                                'Molecular Weight': molecular_weight,
                                'Charge at pH': charge_at_pH,
                                'Aromaticity': aromaticity,
                                'Secondary structure fraction Helix': sec_H,
                                'Secondary structure fraction Turn': sec_T,
                                'Secondary structure fraction Sheet': sec_S,
                                'Boman Index': Boman_index,
                                'Amphipathicity': amphipathicity,
                                'Correlation': correlation,
                                'Covariance': covariance,
                                'Mass': mass,
                                'Mz': mz,
                            }
                            all_features.append(features)

                        except Exception as e:
                            st.error(f"Error parsing Sequence {rec.id}: {e}")

            # -------- TXT --------
            elif file_type.endswith(".txt"):
                file_content = uploaded_file.read().decode("utf-8")

                seq_list = [s.strip() for s in file_content.splitlines() if s.strip()]

                for i, seq in enumerate(seq_list):
                    try:
                        bio_analysis = Bio_analysis(seq)

                        Gravy = bio_analysis.get_gravy()
                        instability_index = bio_analysis.get_instability_index()
                        Aliphatic_Index = bio_analysis.get_aliphatic_index()
                        Boman_index = bio_analysis.get_boman_index()
                        isoelectric_point = bio_analysis.get_isoelectric_point()
                        net_charge = bio_analysis.get_net_charge()
                        molecular_weight = bio_analysis.get_molecular_weight()
                        charge_at_pH = bio_analysis.get_charge_at_pH()
                        aromaticity = bio_analysis.get_aromaticity()
                        sec_H, sec_T, sec_S = bio_analysis.get_secondary_structure_fraction()
                        amphipathicity = bio_analysis.get_amphipathicity()
                        correlation = bio_analysis.get_auto_correlation()
                        covariance = bio_analysis.get_auto_covariance()
                        hydrophobic_moenet = bio_analysis.get_hydrophobic_moenet()
                        mass = bio_analysis.get_mass()
                        mz = bio_analysis.get_mz()
                        SequenceLength = bio_analysis.get_sequenceLength()

                        features = {
                            "Sequence": seq,
                            'Length': SequenceLength,
                            'Gravy': Gravy,
                            'Instability Index': instability_index,
                            'Aliphatic Index': Aliphatic_Index,
                            'Isoelectric point': isoelectric_point,
                            'Net charge': net_charge, 
                            'Molecular Weight': molecular_weight,
                            'Charge at pH': charge_at_pH,
                            'Aromaticity': aromaticity,
                            'Secondary structure fraction Helix': sec_H,
                            'Secondary structure fraction Turn': sec_T,
                            'Secondary structure fraction Sheet': sec_S,
                            'Boman Index': Boman_index,
                            'Amphipathicity': amphipathicity,
                            'Correlation': correlation,
                            'Covariance': covariance,
                            'Mass': mass,
                            'Mz': mz,
                        }
                        all_features.append(features)

                    except Exception as e:
                        st.error(f"Error parsing Sequence line {i+1}: {e}")
        all_features = pd.DataFrame(all_features)

        # ----------------------------
        # Optimization settings
        # ----------------------------
        st.markdown("---")
        st.header("Select Algorithms")
        algorithms = st.multiselect(
            "Choose optimization algorithms",
            options=[
                'NSGA-II', 'NSGA-III', 'R-NSGA-II', 'R-NSGA-III',
                'U-NSGA-III', 'AGE-MOEA', 'AGE-MOEA-II'
            ]
        )

        pop_size = st.number_input("Population size", min_value=10, max_value=400, value=80, step=10)
        length = st.number_input("Peptide Sequence length", min_value=10, max_value=20, value=10, step=1)
        generations = st.number_input("Number of generations", min_value=10, max_value=200, value=100, step=1)
        
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

    if len(Bacteria) == 0 and uploaded_file is None:
        st.warning("Please select at least one optimization algorithm or upload a dataset.")
    
    elif len(opt) < 2:
        st.warning("Please select at least two objectives to optimize.")
    else:
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

        if "constraints" not in st.session_state:
            st.session_state.constraints = []

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

        if st.button("➕ Add Constraint"):
            new_constraint = {
                "Feature": constraint_feature,
                "Type": "max" if "≤" in constraint_type else "min",
                "Value": constraint_value,
            }

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

        if st.button("🗑️ Clear All Constraints"):
            st.session_state.constraints = []
            st.info("All constraints have been cleared.")

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

        if len(Bacteria) == 0 and uploaded_file is None:
            st.warning("Please select at least one bacteria from the sidebar or upload a dataset.")
            st.stop()
        elif uploaded_file is not None:
            st.subheader("Uploaded peptide data preview")
            st.dataframe(all_features)
            df = pd.DataFrame(all_features)
            
            if len(Bacteria) > 0:
                # 合併選擇的細菌 CSV
                dfs = [df]
                for b in Bacteria:
                    try:
                        temp_df = pd.read_csv(f"dataset/biopython-{b}.csv").drop(columns=['MIC'])
                        dfs.append(temp_df)
                    except FileNotFoundError as e:
                        st.error(f"Missing file for {b}: {e}")
                df = pd.concat(dfs, ignore_index=True).drop_duplicates(subset=['Sequence'])
                st.success(f"Merged uploaded data with selected bacteria: {', '.join(Bacteria)}")
            st.dataframe(df)
            
        else:
            # Initialize session state for dataframe
            if "loaded_df" not in st.session_state:
                st.session_state.loaded_df = None
            
            # Load data only if not cached
            if st.session_state.loaded_df is None:
                try:
                    df = pd.read_csv(f"dataset\\biopython-{Bacteria[0]}.csv")
                    for i in range(1, len(Bacteria)):
                        temp_df = pd.read_csv(f"dataset\\biopython-{Bacteria[i]}.csv").drop(columns=['MIC'])
                        df = pd.concat([df, temp_df], ignore_index=True).drop_duplicates(subset=['Sequence'])
                    st.session_state.loaded_df = df
                    st.success(f"Loaded data for {', '.join(Bacteria)} with {len(df)} Sequences.")
                except FileNotFoundError as e:
                    st.error(f"Missing file: {e}")
                    st.stop()
            else:
                df = st.session_state.loaded_df
                st.success(f"Using cached data for {', '.join(Bacteria)} with {len(df)} Sequences.")

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
        else:
            st.stop()
        # ----------------------------
        # Display cached results
        # ----------------------------
        st.subheader("📋 Optimization Results")

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
            if algo not in st.session_state.optimization_results:
                #st.warning(f"Skip {algo} because no results found.")
                continue

            if "merged_df_flipped" not in locals() or merged_df_flipped is None:
                #st.warning(f"Skip {algo} because merged_df_flipped is not defined.")
                continue

            if "pareto_df_flipped" not in locals() or pareto_df_flipped is None:
                #st.warning(f"Skip {algo} because pareto_df_flipped is not defined.")
                continue
            # --------------------------------
            fasta_path = os.path.join(user_home, f"{algo}.fasta")
            with open(fasta_path, "w") as fasta_file:
                for _, row in merged_df_flipped.iterrows():
                    fasta_file.write(f">{row['Sequence']}\n{row['Sequence']}\n")

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
    If you use this app for your research, please cite the following paper:\n
    Yang C-H, Chen Y-L, Cheung T-H, Chuang L-Y. Multi-Objective Optimization Accelerates the De Novo Design of Antimicrobial Peptide for Staphylococcus aureus. International Journal of Molecular Sciences. 2024; 25(24):13688. https://doi.org/10.3390/ijms252413688
    """)

with main_tab3:
    st.header("How to Design Antimicrobial Peptides using this App")
    st.markdown("""
    This app integrates multi-objective optimization frameworks with peptide physicochemical profiling.  
    Follow the workflow below to construct customized antimicrobial peptide (AMP) candidates.

    ### 1. **Select Target Bacteria**
    Choose one or more pathogens from the sidebar.  
    The app will automatically load precomputed physicochemical property datasets associated with the selected species.

    ### 2. **Define Optimization Parameters**
    - **Population Size**: Controls the diversity of candidate sequences within the evolutionary search.
    - **Peptide Length**: Specifies the length of generated sequences for de novo design.
    - **Generations**: Determines the number of iterations for algorithmic evolution.

    These parameters directly influence search convergence, diversity maintenance, and computational runtime.

    ### 3. **Choose Optimization Algorithms**
    Select one or multiple multi-objective evolutionary algorithms (MOEAs), such as:
    - NSGA-II / NSGA-III  
    - R-NSGA-II / R-NSGA-III  
    - U-NSGA-III  
    - AGE-MOEA / AGE-MOEA-II  

    Each algorithm presents different strengths in balancing exploration and convergence toward high-quality Pareto-optimal peptides.

    """)
    st.image("螢幕擷取畫面 2025-11-16 012845.png")
    st.markdown("""
    ### 4. **Select Physicochemical Objectives**
    Choose at least two descriptors for optimization.  
    Available objectives include:
    - Hydrophobicity (Gravy)
    - Instability Index  
    - Isoelectric Point  
    - Net Charge  
    - Aliphatic Index  
    - Aromaticity  
    - Molecular Weight  
    - Boman Index  
    - Secondary Structure Fractions (Helix, Turn, Sheet)

    For each objective, specify whether the algorithm should **minimize** or **maximize** the descriptor.  
    (Hydrophobicity is treated as *hydrophilicity* or *hydrophobicity* depending on user preference.)
                """)
    st.image("螢幕擷取畫面 2025-11-16 013006.png")
    st.markdown("""
    ### 5. **Add Optional Constraints**
    Users may define upper/lower bounds to restrict the peptide search space.  
    For example:
    - Instability Index ≥ 40  
    - 1 ≤ Net Charge
    - Hydrophobicity ≥ 1

    Constraints help enforce biologically realistic design regions and improve hit quality.
                """)
    st.image("螢幕擷取畫面 2025-11-16 013025.png")
    st.markdown("""
    ### 6. **Run Optimization**
    Press **“Run Optimization”** to execute the selected algorithms.  
    The app will:
    - Perform evolutionary optimization  
    - Compute Pareto fronts  
    - Identify non-dominated peptide solutions  
    - Cache results for further analysis  

    ### 7. **View and Download Results**
    Results include:
    - Objective value tables
    - Full Pareto-optimal peptide lists
    - Merged physicochemical property profiles
    - FASTA files for external analysis  
    - CSV exports for downstream modeling

    ### 8. **Visualize Optimization Outcomes**
    You may generate:
    - Multi-dimensional Pareto front plots  
    - Amino acid composition heatmaps  
    - Sequence-level physicochemical distribution analyses  

    These visualizations provide insight into peptide behavior, trade-offs among descriptors, and optimization dynamics.

                """)
    st.image("螢幕擷取畫面 2025-11-16 013229.png")
    st.markdown("""
    ### Summary
    This app provides a structured, multi-objective approach to AMP design by integrating algorithmic search, physicochemical evaluation, and biological constraint modeling. It aims to accelerate the rational development of antimicrobial peptides with optimized properties.
    If you encounter any issues or have questions, please upload your issues or figure to the 
                https://github.com/ksuee108/amp_desige/issues, and we will get back to you as soon as possible.
    """)

with main_tab4:
    st.header("Related Databases and Prediction Websites")
    AMP_databases = {
        "Website":["Peptaibols", "Cybase", "BACTIBASE", "CAMP", "DADP", "HIPdb", "Hemolytik", "ParaPep", "CancerPPD/AntiCP 2.0", "DBAASP", "BaAMPs", "SATPdb", "DRAMP", "InverPep", "ARA-PEPs", "MBPDB", "AntiTbPdb", "LABiocin", "ADAPTABLE", "FoldamerDB", "AntiCP 2.0", "FermFooDb", "B-AMP", "SuPepMem", "ACovPepDB", "AMPDB v1", "DRAVP", "GtoPdb", "aSynPEP-DB", "AbAMPdb", "AVR/I/SSAPDB", "TAMRSA", "IAMPDB", "ABPDB"],
        "Link":["https://peptaibol.cryst.bbk.ac.uk/home.shtml", "https://www.cybase.org.au/", "https://bactibase.pfba-lab-tun.org/main.php", "http://www.bicnirrh.res.in/antimicrobial", "http://split4.pmfst.hr/dadp/", "http://crdd.osdd.net/servers/hipdb/", "http://crdd.osdd.net/raghava/hemolytik/", "http://crdd.osdd.net/raghava/parapep/", "http://crdd.osdd.net/raghava/cancerppd/", "http://dbaasp.org/home.xhtml", "http://www.baamps.it/", "http://crdd.osdd.net/raghava/satpdb/", "http://dramp.cpu-bioinfor.org/", "http://ciencias.medellin.unal.edu.co/gruposdeinvestigacion/prospeccionydisenobiomoleculas/InverPep/public/home_en", "http://www.biw.kuleuven.be/CSB/ARA-PEPs", "http://webs.iiitd.edu.in/raghava/antitbpdb/", "https://mbpdb.nws.oregonstate.edu/", "https://labiocin.univ-lille.fr/", "http://gec.u-picardie.fr/adaptable", "http://foldamerdb.ttk.hu/", "https://webs.iiitd.edu.in/raghava/anticp2/", "https://webs.iiitd.edu.in/raghava/fermfoodb/", "https://b-amp.karishmakaushiklab.com/", "https://supepmem.com/", "http://i.uestc.edu.cn/ACovPepDB/", "https://bblserver.org.in/ampdb/", "http://dravp.cpu-bioinfor.org/", "https://www.guidetopharmacology.org", "https://asynpepdb.ppmclab.com/", "https://abampdb.mgbio.tech/", "https://bblserver.org.in/avrissa/", "https://bblserver.org.in/tamrsar/", "https://bblserver.org.in/iampdb/", "http://www.acdb.plus/ABPDB"],
    }

    st.markdown("### 1. AMP Databases")
    df_AMP_databases = pd.DataFrame(AMP_databases)
    st.table(df_AMP_databases)

    st.markdown("### 2. AMP Prediction Websites")
    AMP_prediction_websites = {
        "Website":["BAGLE", "AntiBP", "AMPer", "CAMP Prediction", "antiSMASH", "AMPA", "AMP_Scanner", "AI4AXP"],
        "Link":["http://bagel.molgenrug.nl/", "https://webs.iiitd.edu.in/raghava/antibp/submit.html", "http://marray.cmdr.ubc.ca/cgi-bin/amp.pl", "http://www.camp.bicnirrh.res.in/predict/", "http://antismash.secondarymetabolites.org/", "http://tcoffee.crg.cat/apps/ampa/do", "http://www.ampscanner.com", "https://axp.iis.sinica.edu.tw/"]
    }
    df_AMP_prediction_websites = pd.DataFrame(AMP_prediction_websites)
    st.table(df_AMP_prediction_websites)