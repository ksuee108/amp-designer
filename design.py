import os
import pandas as pd
import numpy as np

from pymoo.optimize import minimize
from pymoo.core.problem import Problem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.moo.rnsga2 import RNSGA2
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.algorithms.moo.rnsga3 import RNSGA3
from pymoo.algorithms.moo.unsga3 import UNSGA3
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.algorithms.moo.age import AGEMOEA
from pymoo.algorithms.moo.age2 import AGEMOEA2

from sklearn.preprocessing import MinMaxScaler
from pymoo.visualization.pcp import PCP
from pymoo.visualization.petal import Petal
from pymoo.visualization.radviz import Radviz
from pymoo.visualization.scatter import Scatter

from BioAnalysis import Bio_analysis
import matplotlib.pyplot as plt
import streamlit as st

class MyProblemWithData(Problem):
    def __init__(self, sequences, optimization_directions=None, sequence_length=12, opt=2, constraint_dict_list=None, *args, **kwargs):
        super().__init__(n_var=sequence_length, n_obj=len(opt), n_constr=6, xl=0, xu=20)
        self.sequences = sequences
        self.pareto_sequences = []
        self.optimization_directions = optimization_directions
        self.sequence_length = sequence_length
        self.opt = opt
        self.constraint_dict_list = constraint_dict_list
        self.combinations = []

    def _evaluate(self, X, out, *args, **kwargs):
        X = X.astype(np.float64)
        n_samples = X.shape[0]
        objectives = np.zeros((n_samples, len(self.opt)), dtype=np.float64)
        constraints = np.zeros((n_samples, 6), dtype=np.float64)
        aa_str = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
        new_combinations = [''.join([aa_str[int(xi)] for xi in X[i]]) for i in range(n_samples)]
        self.combinations.extend(new_combinations)

        for i, peptide in enumerate(new_combinations):
            bio_analysis = Bio_analysis(peptide)
 
            Gravy = bio_analysis.get_gravy()
            instability_index = bio_analysis.get_instability_index()
            Aliphatic_Index = bio_analysis.get_aliphatic_index()
            Boman_Index = bio_analysis.get_boman_index()
            isoelectric_point = bio_analysis.get_isoelectric_point()
            net_charge = bio_analysis.get_net_charge()
            molecular_weight = bio_analysis.get_molecular_weight()
            charge_at_pH = bio_analysis.get_charge_at_pH()
            aromaticity = bio_analysis.get_aromaticity()
            secondary_structure_fraction_Helix = bio_analysis.get_secondary_structure_fraction()[0]
            secondary_structure_fraction_Turn = bio_analysis.get_secondary_structure_fraction()[1]
            secondary_structure_fraction_Sheet = bio_analysis.get_secondary_structure_fraction()[2]
            sequenceLength = bio_analysis.get_sequenceLength()
            amphipathicity = bio_analysis.get_amphipathicity()
            correlation = bio_analysis.get_auto_correlation()
            covariance = bio_analysis.get_auto_covariance()
            hydrophobic_moenet = bio_analysis.get_hydrophobic_moenet()
            mass = bio_analysis.get_mass()
            mz = bio_analysis.get_mz()
            

            physicochemical_properties = {
            'Gravy': Gravy,
            'Instability Index': instability_index,
            'Aliphatic Index': Aliphatic_Index,
            'Isoelectric point': isoelectric_point,
            'Net charge': net_charge, 
            'Molecular Weight': molecular_weight,
            'Charge at pH': charge_at_pH,
            'Aromaticity': aromaticity,
            'Secondary structure fraction Helix': secondary_structure_fraction_Helix,
            'Secondary structure fraction Turn': secondary_structure_fraction_Turn,
            'Secondary structure fraction Sheet': secondary_structure_fraction_Sheet,
            'Boman Index': Boman_Index,
            'Amphipathicity': amphipathicity,
            'Correlation': correlation,
            'Covariance': covariance,
            'Hydrophobic Moenet': hydrophobic_moenet,
            'Mass': mass,
            'Mz': mz,}
            
            # Set objective values according to selected objectives
            for obj_idx, obj_name in enumerate(self.opt):
                value = physicochemical_properties[obj_name]
                if self.optimization_directions[obj_name] == 'Minimize':
                    objectives[i, obj_idx] = value
                else:  # Maximize
                    objectives[i, obj_idx] = -value

            # Store objective values and sequence
            self.pareto_sequences.append({
                'Sequence': peptide,
                'Gravy': Gravy,
                'Instability Index': instability_index,
                'Aliphatic Index': Aliphatic_Index,
                'Isoelectric point': isoelectric_point,
                'Net charge': net_charge, 
                'Molecular Weight': molecular_weight,
                'Charge at pH': charge_at_pH,
                'Aromaticity': aromaticity,
                'Secondary structure fraction Helix': secondary_structure_fraction_Helix,
                'Secondary structure fraction Turn': secondary_structure_fraction_Turn,
                'Secondary structure fraction Sheet': secondary_structure_fraction_Sheet,
                'Boman Index': Boman_Index,
                'Amphipathicity': amphipathicity,
                'Correlation': correlation,
                'Covariance': covariance,
                'Hydrophobic Moenet': hydrophobic_moenet,
                'Mass': mass,
                'Mz': mz,
                'Sequence Length': sequenceLength
            })
            
            # Constraints setup (example)
            if self.constraint_dict_list:
                for j, constraint in enumerate(self.constraint_dict_list):
                    feature_name = constraint.get("Feature", "Gravy")
                    constraint_type = constraint.get("Type", "max")
                    constraint_value = float(constraint.get("Value"))

                    objective = float(physicochemical_properties.get(feature_name))

                    if constraint_type == "max":
                        constraints[i, j] = constraint_value - objective
                    else:
                        constraints[i, j] = objective - constraint_value
            else:
                pass
        out["F"] = objectives
        out["G"] = constraints

def amino_acid_percentage(path, algorithms):
    amino_acid_counts = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
                         'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

    total_amino_acids = 0

    optimized_result = {}
    
    # Load results for all algorithms
    for name in algorithms:
        filepath = os.path.join(path, f"{name} optimize result.csv")
        optimized_result[name] = pd.read_csv(filepath)
            
    for algo_name, pareto_sequences in optimized_result.items():
        for seq in pareto_sequences['Sequence']:
            for aa in seq:
                if aa in amino_acid_counts:
                    amino_acid_counts[aa] += 1
                    total_amino_acids += 1

        percentages = {aa: count / total_amino_acids * 100 for aa, count in amino_acid_counts.items()}

        amino_acids = list(percentages.keys())
        percentage_values = list(percentages.values())
        fig, ax = plt.subplots()
        ax.bar(amino_acids, percentage_values)
        ax.set_title(f'Percentage of Each Amino Acid {algo_name}')
        ax.set_xlabel('Amino Acid')
        ax.set_ylabel('Percentage')
        fig.tight_layout()
        st.pyplot(fig)
        plt.close(fig)

def plot_pareto_fronts_many(path, algorithms, optimization_directions):
    """
    Plot Pareto fronts visualization for all algorithms.
    
    Args:
        algorithms: list of algorithm names (e.g., ['NSGA-III', 'UNSGA3'])
        optimization_directions: dict of objectives (e.g., {'Gravy': 'Minimize', ...})
        path: directory where result files are stored
    """
    optimized_result = {}
    all_optimized_result = {}
    
    # Load results for all algorithms
    for name in algorithms:
        try:
            filepath = os.path.join(path, f"{name} optimize result.csv")

            if os.path.exists(filepath):
                optimized_result[name] = pd.read_csv(filepath)
            else:
                st.warning(f"⚠️ Main optimization result not found for {name}: {filepath}")
                continue

            filepath_all_result = os.path.join(path, f"{name} all optimize result.csv")
            if os.path.exists(filepath_all_result):
                all_optimized_result[name] = pd.read_csv(filepath_all_result)
            else:
                st.info(f"ℹ️ No 'all results' file found for {name} (skipped).")

        except pd.errors.EmptyDataError:
            st.error(f"❌ CSV file for {name} is empty or corrupted: {filepath}")
            continue

        except Exception as e:
            st.error(f"❌ Error loading results for {name}: {e}")
            continue

    selected_cols_name = list(optimization_directions.keys())
    st.write(f"**Objectives:** {selected_cols_name}")
    
    for algo_name, df in optimized_result.items():
        all_result_df = all_optimized_result[algo_name]
        if 'Gravy' in selected_cols_name and optimization_directions['Gravy']=="hydrophobicity":
            all_result_df['Gravy'] = -all_optimized_result[algo_name]['Gravy']

        for feature, mode in optimization_directions.items():
            if optimization_directions[feature]=="Maximize":
                all_result_df[feature] = -all_optimized_result[algo_name][feature]
        st.write(f"### {algo_name}")
        
        missing_cols = [col for col in selected_cols_name if col not in df.columns]
        if missing_cols:
            st.warning(f"Missing columns in {algo_name}: {missing_cols}")
            continue
        
        X = df[selected_cols_name].values
        n_samples = X.shape[0]
        n_objectives = len(selected_cols_name)

        all_X = all_result_df[selected_cols_name].values
        
        if n_samples == 0:
            st.warning(f"No data in {algo_name}")
            continue
        
        # Normalize data
        scaler = MinMaxScaler()
        X_scaled = scaler.fit_transform(X)
        all_X_scaled = scaler.fit_transform(all_X)
        
        # Color map
        colors = plt.cm.viridis(np.linspace(0, 1, n_samples))
        colors = [tuple(c) for c in colors]
        
        try:
            # Parallel coordinates plot (matplotlib)
            st.write("**Parallel Coordinates Plot:**")
            fig, ax = plt.subplots(figsize=(12, 6))
            
            # Draw each line
            for idx in range(n_samples):
                ax.plot(range(len(selected_cols_name)), X_scaled[idx, :], 
                       color=colors[idx], alpha=0.6, linewidth=1)
            
            ax.set_xticks(range(len(selected_cols_name)))
            ax.set_xticklabels(selected_cols_name, rotation=45, ha='right')
            ax.set_ylabel('Normalized Value')
            ax.set_title(f'Parallel Coordinates - {algo_name}')
            ax.grid(True, alpha=0.3)
            st.pyplot(fig)
            plt.close(fig)
            
            # Scatter matrix plot
            st.write("**Scatter Matrix Plot:**")
            fig, axes = plt.subplots(n_objectives, n_objectives, figsize=(14, 14))
            
            # If single objective, make axes 2D array
            if n_objectives == 1:
                axes = np.array([[axes]])
            elif n_objectives > 1:
                if axes.ndim == 1:
                    axes = axes.reshape(-1, 1)
            with st.spinner("Running optimization... This may take a few minutes."):
                for i in range(n_objectives):
                    for j in range(n_objectives):
                        ax = axes[i, j] if n_objectives > 1 else axes[0, 0]
                        
                        if i == j:
                            # Show objective name on diagonal
                            ax.text(0.5, 0.5, selected_cols_name[i], 
                                ha='center', va='center', fontsize=12, fontweight='bold')
                            ax.set_xticks([])
                            ax.set_yticks([])
                            ax.spines['top'].set_visible(False)
                            ax.spines['right'].set_visible(False)
                            ax.spines['bottom'].set_visible(False)
                            ax.spines['left'].set_visible(False)
                        else:
                            # Scatter plot for non-diagonal
                            scatter = ax.scatter(all_X[:, j], all_X[:, i], 
                                            c=np.arange(len(all_X)), cmap='viridis', 
                                            s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
                            ax.scatter(X[:, j], X[:, i], 
                                            color='red', s=50, alpha=0.7, edgecolors='red', linewidth=0.5)
                            
                            # Labels
                            if j == 0:
                                ax.set_ylabel(selected_cols_name[i], fontsize=10, fontweight='bold')
                            else:
                                ax.set_ylabel('')
                            
                            if i == n_objectives - 1:
                                ax.set_xlabel(selected_cols_name[j], fontsize=10, fontweight='bold')
                            else:
                                ax.set_xlabel('')
                            ax.grid(True, alpha=0.3)
                
            fig.suptitle(f'Scatter Matrix Plot - {algo_name}', fontsize=16, fontweight='bold')
            plt.tight_layout()
            st.pyplot(fig)
            plt.close(fig)
        
        except Exception as e:
            st.error(f"Error plotting {algo_name}: {e}")
            continue

def plot_pareto_fronts_multi(path, algorithms, optimization_directions):
    """
    Plot Pareto fronts visualization for all algorithms.
    This function supports 2D, 3D and scatter matrix (for 3+ objectives).
    """
    optimized_result = {}
    all_optimized_result = {}
    
    # Load results for all algorithms
    for name in algorithms:
        try:
            filepath = os.path.join(path, f"{name} optimize result.csv")
            if os.path.exists(filepath):
                optimized_result[name] = pd.read_csv(filepath)
            else:
                st.warning(f"⚠️ Main optimization result not found for {name}: {filepath}")
                continue

            filepath_all_result = os.path.join(path, f"{name} all optimize result.csv")
            if os.path.exists(filepath_all_result):
                all_optimized_result[name] = pd.read_csv(filepath_all_result)
            else:
                st.info(f"ℹ️ No 'all results' file found for {name} (skipped).")

        except pd.errors.EmptyDataError:
            st.error(f"❌ CSV file for {name} is empty or corrupted: {filepath}")
            continue

        except Exception as e:
            st.error(f"❌ Error loading results for {name}: {e}")
            continue

    selected_cols_name = list(optimization_directions.keys())
    st.write(f"**Objectives:** {selected_cols_name}")
    
    for algo_name, df in optimized_result.items():
        all_result_df = all_optimized_result[algo_name]
        if 'Gravy' in selected_cols_name and optimization_directions['Gravy']=="hydrophobicity":
            all_result_df['Gravy'] = -all_optimized_result[algo_name]['Gravy']

        for feature, mode in optimization_directions.items():
            if optimization_directions[feature]=="Maximize":
                all_result_df[feature] = -all_optimized_result[algo_name][feature]
        st.write(f"### {algo_name}")
        
        missing_cols = [col for col in selected_cols_name if col not in df.columns]
        if missing_cols:
            st.warning(f"Missing columns in {algo_name}: {missing_cols}")
            continue
        
        X = df[selected_cols_name].values
        n_samples = X.shape[0]
        n_objectives = len(selected_cols_name)

        all_X = all_result_df[selected_cols_name].values
        
        if n_samples == 0:
            st.warning(f"No data in {algo_name}")
            continue
        
        try:
            # 2D plot
            if n_objectives == 2:
                st.write("**2D Scatter Plot:**")
                fig, ax = plt.subplots(figsize=(8, 6))
                scatter = ax.scatter(all_X[:, 0], all_X[:, 1], c=np.arange(len(all_X)), 
                                    cmap='viridis', s=10, alpha=0.7, edgecolors='black', linewidth=0.5)
                ax.scatter(X[:, 0], X[:, 1], color='red', s=10, alpha=0.7, edgecolors='red', linewidth=0.5)
                
                ax.set_xlabel(selected_cols_name[0], fontsize=12, fontweight='bold')
                ax.set_ylabel(selected_cols_name[1], fontsize=12, fontweight='bold')
                ax.set_title(f'2D Pareto Front - {algo_name}', fontsize=14, fontweight='bold')
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label('Sample Index', fontsize=11)
                ax.grid(True, alpha=0.3)
                st.pyplot(fig)
                plt.close(fig)
            
            # 3D plot
            elif n_objectives == 3:
                st.write("**3D Scatter Plot:**")
                from mpl_toolkits.mplot3d import Axes3D
                fig = plt.figure(figsize=(10, 8))
                ax = fig.add_subplot(111, projection='3d')
                scatter = ax.scatter(all_X[:, 0], all_X[:, 1], all_X[:, 2], 
                                   c=np.arange(len(all_X)), cmap='viridis', s=10, alpha=0.7, edgecolors='black', linewidth=0.5)
                ax.scatter(X[:, 0], X[:, 1], X[:, 2], 
                                   color='red', s=10, alpha=0.7, edgecolors='red', linewidth=0.5)
                ax.set_xlabel(selected_cols_name[0], fontsize=10, fontweight='bold')
                ax.set_ylabel(selected_cols_name[1], fontsize=10, fontweight='bold')
                ax.set_zlabel(selected_cols_name[2], fontsize=10, fontweight='bold')
                ax.set_title(f'3D Pareto Front - {algo_name}', fontsize=12, fontweight='bold')
                cbar = plt.colorbar(scatter, ax=ax, shrink=0.5)
                cbar.set_label('Sample Index', fontsize=10)
                st.pyplot(fig)
                plt.close(fig)
            
            # For 3 or more objectives show scatter matrix
            if n_objectives >= 3:
                st.write("**Scatter Matrix Plot:**")
                
                scaler = MinMaxScaler()
                X_scaled = scaler.fit_transform(X)
                
                fig, axes = plt.subplots(n_objectives, n_objectives, figsize=(14, 14))
                
                if n_objectives == 1:
                    axes = np.array([[axes]])
                elif n_objectives > 1:
                    if axes.ndim == 1:
                        axes = axes.reshape(-1, 1)
                
                for i in range(n_objectives):
                    for j in range(n_objectives):
                        ax = axes[i, j] if n_objectives > 1 else axes[0, 0]
                        
                        if i == j:
                            ax.text(0.5, 0.5, selected_cols_name[i], 
                                ha='center', va='center', fontsize=12, fontweight='bold')
                            ax.set_xticks([])
                            ax.set_yticks([])
                            ax.spines['top'].set_visible(False)
                            ax.spines['right'].set_visible(False)
                            ax.spines['bottom'].set_visible(False)
                            ax.spines['left'].set_visible(False)
                        else:
                            scatter = ax.scatter(X_scaled[:, j], X_scaled[:, i], 
                                            c=np.arange(n_samples), cmap='viridis', 
                                            s=50, alpha=0.7, edgecolors='black', linewidth=0.5)
                            
                            if j == 0:
                                ax.set_ylabel(selected_cols_name[i], fontsize=10, fontweight='bold')
                            else:
                                ax.set_ylabel('')
                            
                            if i == n_objectives - 1:
                                ax.set_xlabel(selected_cols_name[j], fontsize=10, fontweight='bold')
                            else:
                                ax.set_xlabel('')
                            
                            ax.set_xlim([-0.05, 1.05])
                            ax.set_ylim([-0.05, 1.05])
                            ax.grid(True, alpha=0.3)
                
                fig.suptitle(f'Scatter Matrix Plot - {algo_name}', fontsize=16, fontweight='bold')
                plt.tight_layout()
                st.pyplot(fig)
                plt.close(fig)
            
        except Exception as e:
            st.error(f"Error plotting {algo_name}: {e}")
            print(f"Error details: {e}")
            continue

# ----------------------------
class algorithms_setup():
    def __init__(self, path, df, algorithms_list, pop_size, generations, optimization_directions, length, opt, constraint_dict_list):
        self.path = path
        self.sequences = df['Sequence']
        self.algorithms_list = algorithms_list
        self.pop_size = pop_size
        self.generations = generations
        self.optimization_directions = optimization_directions
        self.length = length
        self.opt = opt
        self.constraint_dict_list = constraint_dict_list

    def run_optimization(self):
        # Initialize problem object
        self.problem = MyProblemWithData(
                sequences=self.sequences,
                optimization_directions=self.optimization_directions,
                sequence_length=self.length,
                opt=self.opt,
                constraint_dict_list=self.constraint_dict_list
            )
        st.write("The package program is requesting that multiply an objective that is supposed to be maximized by −1 and minimized for valuation.")
        self.pf = self.problem.pareto_front()

        # reference directions
        ref_dirs = get_reference_directions("das-dennis", n_dim=len(self.opt), n_partitions = 50)

        n_ref_points = len(ref_dirs)

        # Suggested pop_per_ref_point
        suggested_pop_per_ref = 500 // n_ref_points
        if suggested_pop_per_ref < 1:
            suggested_pop_per_ref = 1

        try:
            
            for feature, mode in self.optimization_directions.key():
                if feature == "Gravy":
                    if mode.lower() in ["hydrophobicity"]:
                        self.problem.flip_objective(feature)
                else:
                    if mode.lower() == "maximize":
                        self.problem.flip_objective(feature)
                        
            
            all_algorithms = [[NSGA2(pop_size=self.pop_size),'NSGA-II'],
                              [NSGA3(ref_dirs = ref_dirs,pop_size=self.pop_size),'NSGA-III'], 
                              [UNSGA3(ref_dirs = ref_dirs,pop_size=self.pop_size),'U-NSGA-III'], 
                              [RNSGA2(ref_points = ref_dirs, pop_size=self.pop_size),'R-NSGA-II'], 
                              [RNSGA3(ref_points=ref_dirs, pop_per_ref_point=self.pop_size), 'R-NSGA-III'], 
                              [AGEMOEA(pop_size=self.pop_size), 'AGE-MOEA'], 
                              [AGEMOEA2(pop_size=self.pop_size), 'AGE-MOEA-II']]
            self.algorithms = [algo for algo in all_algorithms if algo[1] in self.algorithms_list]
            print(self.algorithms)
        except:
            all_algorithms = [[NSGA2(pop_size=self.pop_size),'NSGA-II'],
                               [NSGA3(ref_dirs = ref_dirs,pop_size=self.pop_size),'NSGA-III'], 
                               [UNSGA3(ref_dirs = ref_dirs,pop_size=self.pop_size),'U-NSGA-III'], 
                               [RNSGA2(ref_points = ref_dirs, pop_size=self.pop_size),'R-NSGA-II'], 
                               [RNSGA3(ref_points=ref_dirs, pop_per_ref_point=suggested_pop_per_ref), 'R-NSGA-III'],
                               [AGEMOEA(pop_size=self.pop_size), 'AGE-MOEA'], 
                               [AGEMOEA2(pop_size=self.pop_size), 'AGE-MOEA-II']]
            self.algorithms = [algo for algo in all_algorithms if algo[1] in self.algorithms_list]
            print(self.algorithms)

    def run(self):
        st.info("Starting optimization...")
        
        # Initialize session state for results
        if "optimization_results" not in st.session_state:
            st.session_state.optimization_results = {}
        
        for algorithm, name in self.algorithms:
            print(name)
            if name in 'NSGA3':
                res = minimize(self.problem, algorithm, seed=1, pf=self.pf, termination=('n_gen', self.generations), verbose=False)
            else:
                res = minimize(self.problem, algorithm, ('n_gen', self.generations), pf=self.pf, verbose=False)
            
            # debug
            # res_dict = pd.DataFrame(res.F)
            # st.dataframe(res_dict)

            res_dict = {}
            for i, obj_name in enumerate(self.optimization_directions.keys()):
                res_dict[obj_name] = res.F[:, i]
            
            res_dict = pd.DataFrame(res_dict)

            pareto_dict = {
                'Sequence': [entry['Sequence'] for entry in self.problem.pareto_sequences],
                'Gravy': [entry['Gravy'] for entry in self.problem.pareto_sequences],
                'Instability Index': [entry['Instability Index'] for entry in self.problem.pareto_sequences],
                'Aliphatic Index': [entry['Aliphatic Index'] for entry in self.problem.pareto_sequences],
                'Boman Index': [entry['Boman Index'] for entry in self.problem.pareto_sequences],
                'Molecular Weight': [entry['Molecular Weight'] for entry in self.problem.pareto_sequences],
                'Sequence Length': [entry['Sequence Length'] for entry in self.problem.pareto_sequences],
                'Isoelectric point': [entry['Isoelectric point'] for entry in self.problem.pareto_sequences],
                'Net charge': [entry['Net charge'] for entry in self.problem.pareto_sequences],
                'Charge at pH': [entry['Charge at pH'] for entry in self.problem.pareto_sequences],
                'Aromaticity': [entry['Aromaticity'] for entry in self.problem.pareto_sequences],
                'Secondary structure fraction Helix': [entry['Secondary structure fraction Helix'] for entry in self.problem.pareto_sequences],
                'Secondary structure fraction Turn': [entry['Secondary structure fraction Turn'] for entry in self.problem.pareto_sequences],
                'Secondary structure fraction Sheet': [entry['Secondary structure fraction Sheet'] for entry in self.problem.pareto_sequences],
                'Amphipathicity': [entry['Amphipathicity'] for entry in self.problem.pareto_sequences],
                'Correlation': [entry['Correlation'] for entry in self.problem.pareto_sequences],
                'Covariance': [entry['Covariance'] for entry in self.problem.pareto_sequences],
                'Hydrophobic Moenet': [entry['Hydrophobic Moenet'] for entry in self.problem.pareto_sequences],
                'Mass': [entry['Mass'] for entry in self.problem.pareto_sequences],
                'Mz': [entry['Mz'] for entry in self.problem.pareto_sequences],
                'model': [name for entry in self.problem.pareto_sequences]
            }

            pareto_df = pd.DataFrame(pareto_dict)
            merged_df = pd.concat([res_dict, pareto_df], axis=1)
            merged_df = merged_df.loc[:,~merged_df.columns.duplicated()].dropna()
            if "optimization_results" not in st.session_state:
                st.session_state.optimization_results = {}

            st.session_state.optimization_results[name] = {
                "res_dict": res_dict,
                "pareto_df": pareto_df,
                "merged_df": merged_df
            }
        return list(st.session_state.optimization_results.keys())