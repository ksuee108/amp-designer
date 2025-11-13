import os
import pandas as pd
import numpy as np

from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
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

from BioAnalysis import Bio_analysis
import matplotlib.pyplot as plt
import tensorflow as tf
#import helixvis
import tensorflow as tf
from keras import mixed_precision
import math
from collections import Counter

mixed_precision.set_global_policy('mixed_float16')
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)


def plot(name):
    plot = Scatter(plot_3d=False, tight_layout=True, labels = ['Gravy', 'Instability Index', 'Aliphatic Index', "Boman_Index"], title = f"{name}", legend = True)
    plot.add(res.F, label="Pareto Front", s=10)
    plot.show()

def amino_acid_percentage():
    amino_acid_counts = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0,
                         'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

    total_amino_acids = 0

    for seq in pareto_df['sequence']:
        for aa in seq:
            if aa in amino_acid_counts:
                amino_acid_counts[aa] += 1
                total_amino_acids += 1

    percentages = {aa: count / total_amino_acids * 100 for aa, count in amino_acid_counts.items()}

    amino_acids = list(percentages.keys())
    percentage_values = list(percentages.values())

    plt.bar(amino_acids, percentage_values)
    plt.title('Percentage of Each Amino Acid')
    plt.xlabel('Amino Acid')
    plt.ylabel('Percentage')
    plt.tight_layout()
    plt.show()

class MyProblemWithData(Problem):
    def __init__(self, combinations):
        super().__init__(n_var=12, n_obj=4, n_constr=6, xl=0, xu=20)
        self.combinations = combinations
        self.pareto_sequences = []

    def _evaluate(self, X, out, *args, **kwargs):
        X = X.astype(np.float64)
        n_samples = X.shape[0]
        objectives = np.zeros((n_samples, 4), dtype=np.float64)
        constraints = np.zeros((n_samples, 6), dtype=np.float64)
        aa_str = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
        
        new_combinations = [''.join([aa_str[int(xi)] for xi in X[i]]) for i in range(n_samples)]
        self.combinations.extend(new_combinations)

        for i, peptide in enumerate(new_combinations):
            bio_analysis = Bio_analysis(peptide)
 
            Gravy = bio_analysis.get_gravy() # 疏水性
            instability_index = bio_analysis.get_instability_index()
            Aliphatic_Index = bio_analysis.get_aliphatic_index()
            Boman_index = bio_analysis.get_boman_index() * -1 # 鮑曼指數
            #net_charge = bio_analysis.get_net_charge()  # 

            # 特徵
            # 计算目标值
            objectives[i, 0] = float(Gravy)
            objectives[i, 1] = float(instability_index)
            objectives[i, 2] = float(Aliphatic_Index)
            objectives[i, 3] = float(Boman_index)


            # 存储目标值和序列
            self.pareto_sequences.append({
                'sequence': peptide,
                'Gravy': Gravy,
                'Instability Index': instability_index,
                'Aliphatic Index': Aliphatic_Index,
                'isoelectric_point': bio_analysis.get_isoelectric_point(),
                'net_charge': bio_analysis.get_net_charge(), 
                'Molecular Weight': bio_analysis.get_molecular_weight(),
                'charge_at_pH': bio_analysis.get_charge_at_pH(),
                'aromaticity': bio_analysis.get_aromaticity(),
                'secondary_structure_fraction_Helix': bio_analysis.get_secondary_structure_fraction()[0],
                'secondary_structure_fraction_Turn': bio_analysis.get_secondary_structure_fraction()[1],
                'secondary_structure_fraction_Sheet': bio_analysis.get_secondary_structure_fraction()[2],
                'Boman_Index': Boman_index,
                'sequenceLength': bio_analysis.get_sequenceLength()
                
            })
            # G = x -bio_analysis... # 大於 x
            # G = bio_analysis... - x # 小於 x
            # 定义约束条件
            constraints[i, 0] = float(1 - bio_analysis.get_gravy())# 疏水性约束
            #constraints[i, 1] = float(net_charge - (-4)) # 淨電荷约束
            #constraints[i, 1] = float(bio_analysis.get_instability_index() - 40)# 不稳定性约束
            #constraints[i, 2] = float(71 - bio_analysis.get_aliphatic_index())# 脂肪族约束
            #constraints[i, 2] = float(bio_analysis.get_boman_index() - 2.5)# 鮑曼约束
            #constraints[i, 4] = float(bio_analysis.get_molecular_weight()  - 600)# 分子量约束
            
        out["F"] = objectives
        out["G"] = constraints

    def save_sequences_to_csv(self, filename):
        df = pd.DataFrame(self.pareto_sequences).drop_duplicates('sequence')
        df.to_csv(os.path.join(path, filename), index=False)

def plot_common_sequences_in_chunks(common_df, selected_cols, chunk_size=9):
    num_sequences = len(common_df)
    print(num_sequences)
    num_chunks = math.ceil(num_sequences / chunk_size)

    for i in range(num_chunks):
        start = i * chunk_size
        end = min((i + 1) * chunk_size, num_sequences)

        chunk = common_df.iloc[start:end]
        titles = [f"{s}" for s in chunk["sequence"]]

        plot = Petal(bounds=[0, 1], title=titles)
        plot.add(common_df[selected_cols].values)
        plot.show()

def plot_pareto_fronts(path):
    unsga3_df = pd.read_csv(os.path.join(path, 'UNSGA3 merged_result1.csv'))
    nsga3_df = pd.read_csv(os.path.join(path, 'NSGA3 merged_result1.csv'))
    #rnsga3_df = pd.read_csv(os.path.join(path, 'RNSGA3 merged_result1.csv'))
    agemoea_df = pd.read_csv(os.path.join(path, 'AGEMOEA merged_result1.csv'))
    agemoea2_df = pd.read_csv(os.path.join(path, 'AGEMOEA2 merged_result1.csv'))
    
    selected_cols = ["Gravy", "Instability Index", "Aliphatic Index", "Boman_Index", "Molecular Weight", "isoelectric_point", "net_charge"]
    selected_cols_name = ["Gravy", "Instability Index", "Aliphatic Index", "Boman Index", "Molecular Weight", "isoelectric point", "net charge"]
    for model in [unsga3_df, nsga3_df, agemoea_df, agemoea2_df]:
        
        X = model[selected_cols].values
        n_samples = X.shape[0]

        scaler = MinMaxScaler()
        X_scaled = scaler.fit_transform(X)

        # 平行座標圖
        colors = plt.cm.viridis(np.linspace(0, 1, n_samples))
        colors = [tuple(c) for c in colors]   # 轉成 [(r,g,b,a), ...]

        # 每條線一個顏色
        pcp = PCP(labels=selected_cols_name)
        rad = Radviz(labels=selected_cols_name)
        for idx, seq in enumerate(model['sequence']):
            pcp.add(X_scaled[idx:idx+1, :], color=colors[idx])   # 每次只加一條線，指定顏色
            rad.add(X_scaled[idx:idx+1, :], color=colors[idx])
        pcp.show()
        rad.show()
    
    selected_cols = ["Gravy", "Instability Index", "Aliphatic Index", "Boman_Index"]
    selected_cols_name = ["Gravy", "Instability Index", "Aliphatic Index", "Boman Index"]

    unsga3_data = unsga3_df[selected_cols].values
    nsga3_data = nsga3_df[selected_cols].values
    #rnsga3_data = rnsga3_df[selected_cols].values
    agemoea_data = agemoea_df[selected_cols].values
    agemoea2_data = agemoea2_df[selected_cols].values
    
    plot2 = Scatter(plot_3d=True, labels = selected_cols, 
                    title = "Pareto Fronts", legend = True)
    #plot2.add(res.F, label="Pareto Front")
    plot2.add(unsga3_data, label="UNSGA3", s=15, alpha=0.8, marker='o')
    plot2.add(nsga3_data, label="NSGA3", s=15, alpha=0.8, marker='^')
    #plot2.add(rnsga3_data, label="RNSGA3", s=15, alpha=0.8, marker='x')
    plot2.add(agemoea_data, label="AGEMOEA", s=15, alpha=0.8, marker='s')
    plot2.add(agemoea2_data, label="AGEMOEA2", s=15, alpha=0.8, marker='p')
    plot2.show()

    common_df = pd.read_csv(os.path.join(path, 'common_sequences1.csv'))
    # 過濾 nsga2_df
    #common_df = nsga3_df[nsga3_df['sequence'].isin(common_sequences)].copy()

    # 翻轉 Gravy 和 Boman_Index
    common_df[['Boman_Index']] *= -1

    # 輸出結果
    common_df.to_csv(os.path.join(path, 'common_sequences2.csv'), index=False)
    
    X = common_df[selected_cols].values
    n_samples = X.shape[0]

    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)

    # 平行座標圖
    colors = plt.cm.viridis(np.linspace(0, 1, n_samples))
    colors = [tuple(c) for c in colors]   # 轉成 [(r,g,b,a), ...]

    # 每條線一個顏色
    pcp = PCP(labels=selected_cols_name)
    rad = Radviz(labels=selected_cols_name)
    for idx, seq in enumerate(common_df['sequence']):
        pcp.add(X_scaled[idx:idx+1, :], color=colors[idx])   # 每次只加一條線，指定顏色
        rad.add(X_scaled[idx:idx+1, :], color=colors[idx])
        
    pcp.show()
    rad.show()

    for idx, seq in enumerate(common_df['sequence']):
        titles = [f"{s}" for s in [seq]]
        plot = Petal(bounds=[0, 1], title=titles)
        plot.add(X_scaled[idx:idx+1, :])
        plot.show()

    # 儲存成 FASTA
    with open(os.path.join(path, 'common_sequences1.fasta'), 'w') as fasta_file:
        for idx, row in common_df.iterrows():
            fasta_file.write(f">{row['sequence']}\n{row['sequence']}\n")

Bacteria =['Escherichia coli', 'Staphylococcus aureus']

file_name_E = 'Biopython-Escherichia coli.csv'
path_E = 'D:\\張\\科技部\\Escherichia coli'
file_name_S = 'Biopython-Staphylococcus aureus.csv'
path_S = 'D:\\張\\科技部\\Staphylococcus aureus'
path = 'D:\\張\\科技部\\com'

data_E = pd.read_csv(os.path.join(path_E, file_name_E))
data_S = pd.read_csv(os.path.join(path_S, file_name_S))
data = pd.concat([data_E, data_S], ignore_index=True)#, data_S
sequences = data['sequence'].tolist()

problem = MyProblemWithData(sequences)
pf = problem.pareto_front()
print(pf)

#參考點
ref_points = np.array([[-2, 2.183333333, 96.66666667, 2.265], [-0.983333333, -84.21666667, 130, 0.07], 
                       [-0.216666667, -73.31666667, 81.66666667, -0.446666667], [-2, -34.11666667, 195, 2.345],
                       [-1.9, -62.41666667, 161.6666667, 2.051666667],[-1.666666667, -48.26666667, 113.3333333, 1.875],
                       [-1.033333333, -73.31666667, 130, 1.001666667], [-0.216666667, -84.21666667, 81.66666667, -0.446666667]])
ref_dirs = get_reference_directions("das-dennis", 4, n_partitions=10, seed=1)  # 使用 das-dennis 方法生成参考方向

n_ref_points = len(ref_dirs)
print(f"參考點數量: {n_ref_points}")

# 計算建議的 pop_per_ref_point
suggested_pop_per_ref = 500 // n_ref_points
if suggested_pop_per_ref < 1:
    suggested_pop_per_ref = 1

print(f"建議的 pop_per_ref_point: {suggested_pop_per_ref}")

try:
    algorithms = [[NSGA2(pop_size=200),'NSGA2'][NSGA3(ref_dirs = ref_dirs,pop_size=500),'NSGA3'], [UNSGA3(ref_dirs = ref_dirs,pop_size=500),'UNSGA3'], 
                  [RNSGA2(ref_points = ref_dirs, pop_size=500),'RNSGA2'], [RNSGA3(ref_points=ref_dirs, pop_per_ref_point=suggested_pop_per_ref), 'RNSGA3'], 
                  [AGEMOEA(pop_size=500), 'AGEMOEA'], [AGEMOEA2(pop_size=500), 'AGEMOEA2']] 
              
except:    
    algorithms = [[NSGA3(ref_dirs = ref_dirs,pop_size=500),'NSGA3'], [AGEMOEA(pop_size=500), 'AGEMOEA'], [AGEMOEA2(pop_size=500), 'AGEMOEA2'], [UNSGA3(ref_dirs = ref_dirs,pop_size=500),'UNSGA3']]
              #,[RNSGA2(ref_points = ref_dirs, pop_size=500),'RNSGA2']]#, , 
              #[NSGA2(pop_size=200),'NSGA2'],

for algorithm, name in algorithms:
    pareto_dict = {}
    print(name)
    if name in 'NSGA3':
        res = minimize(problem, algorithm, seed = 1, pf=pf, termination = ('n_gen', 200), verbose=False)
    else:
        res = minimize(problem, algorithm, ('n_gen', 200), pf=pf, verbose=False)

    res_df = pd.DataFrame(res.F, columns=["Gravy", "Instability Index", "Aliphatic Index", "Boman_Index"])
    res_df.to_csv(os.path.join(path, f"{name} result.csv"), index=False)

    pareto_dict = {
        'sequence': [entry['sequence'] for entry in problem.pareto_sequences],
        'Gravy': [entry['Gravy'] for entry in problem.pareto_sequences],
        'Instability Index': [entry['Instability Index'] for entry in problem.pareto_sequences],
        'Aliphatic Index': [entry['Aliphatic Index'] for entry in problem.pareto_sequences],
        'Boman_Index': [entry['Boman_Index'] for entry in problem.pareto_sequences],
        'Molecular Weight': [entry['Molecular Weight'] for entry in problem.pareto_sequences],
        'sequenceLength': [entry['sequenceLength'] for entry in problem.pareto_sequences],
        'isoelectric_point': [entry['isoelectric_point'] for entry in problem.pareto_sequences],
        'net_charge': [entry['net_charge'] for entry in problem.pareto_sequences],
        'model': [name for entry in problem.pareto_sequences]
    }

    pareto_df = pd.DataFrame(pareto_dict)

    pareto_df.to_csv(os.path.join(path, f"{name} result1.csv"), index=False)

    merged_df = res_df.merge(pareto_df, on=["Gravy", "Instability Index", "Aliphatic Index", "Boman_Index"]).drop_duplicates('sequence')
    merged_df.to_csv(os.path.join(path, f"{name} merged_result1.csv"), index=False)

    # 另存FASTA格式
    with open(os.path.join(path, f"{name}.fasta"), "w") as fasta_file:
        for i, row in merged_df.iterrows():
            fasta_file.write(f">{row['sequence']}\n")
            fasta_file.write(f"{row['sequence']}\n")

plot_pareto_fronts(path)
