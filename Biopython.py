import os ,pandas as pd, numpy as np
from sklearn.model_selection import train_test_split

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.ExPASy import ScanProsite
from Bio import ExPASy
from keras import mixed_precision
import peptides

Bacteria =['Staphylococcus aureus'] #'Pseudomonas aeruginosa', 'Escherichia coli', , 'Acinetobacter baumannii'

for bacteria in Bacteria:
    print(f"處理細菌: {bacteria}")
    path = f'D:\\張\\科技部\\{bacteria}\\bio-{bacteria}.csv'
    data = pd.read_csv(path)
    X = data.drop(columns=['MIC'])#, 'activity', 'unit'])
    y = data['MIC']
    # Split into training and validation sets
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
    sequence_column = 'sequence'
    """
    # 計算序列長度
    X_train_lengths = X_train[sequence_column].apply(len)
    X_val_lengths = X_val[sequence_column].apply(len)

    # 統計描述
    print("X_train 序列長度統計:")
    print(X_train_lengths.describe())
    print("\nX_val 序列長度統計:")
    print(X_val_lengths.describe())

    # 繪製序列長度分佈的直方圖
    plt.figure(figsize=(12, 6))
    plt.hist(data[sequence_column].apply(len), bins=30, color='skyblue', alpha=0.7, label='Train set')
    #plt.hist(X_val_lengths, bins=30, color='lightcoral', alpha=0.7, label='Test set')
    plt.title("Length Distribution")
    plt.xlabel("Length")
    plt.ylabel("Occurrence")
    plt.tight_layout()
    plt.show()"""
    sequenceLength, molecular, isoelectric, charge_pH, Aromaticity, gravy, instability, Helix, Turn, Sheet, net_charge, Aliphatic_Index ,Boman_index,amphipathicity = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
    correlation, covariance, hydrophobic_moment, mass, membrane_position_profile, mz= [],[],[],[],[],[]
    table = peptides.tables.HYDROPHOBICITY["KyteDoolittle"]

    def getNetCharge(seq):
        return sum(1.0 if aa in 'RKH' else -1.0 if aa in 'DE' else 0.0 for aa in seq)

    def getAliphaticIndex(seq):

        A_count, I_count, L_count, V_count = seq.count('A'), seq.count('I'), seq.count('L'), seq.count('V')
        total_length = len(seq)

        A_percent = (A_count / total_length) * 100
        V_percent = (V_count / total_length) * 100
        IL_percent = ((I_count + L_count) / total_length) * 100
        a, b = 2.9, 3.9
        return A_percent + (a * V_percent) + (b * IL_percent)

    def getbomanindex(seq):
        boman = {
            "A": 1.81, "C": 1.28, "D": -8.72, "E": -6.81, "F": 2.98, "G": 0.94,
            "H": -4.66, "I": 4.92, "K": -5.55, "L": 4.92, "M": 2.35, "N": -6.64,
            "P": 0, "Q": -5.54, "R": -14.92, "S": -3.4, "T": -2.57, "V": 4.04,
            "W": 2.33, "Y": -0.14
        }
        return -sum(boman[aa] for aa in seq) / len(seq)

    def calculate_amphipathicity(seq, hydrophobicity_values, hydrophilicity_values):
        hydrophobicity_sum = sum(hydrophobicity_values[aa] for aa in seq)
        hydrophilicity_sum = sum(hydrophilicity_values[aa] for aa in seq)
        return (hydrophobicity_sum / len(seq)) - (hydrophilicity_sum / len(seq))

    for x in data['sequence']:
        print(f"處理細菌: {bacteria}")
        prot_param = ProteinAnalysis(x)
        
        # 计算氨基酸组成
        aa_composition = prot_param.get_amino_acids_percent()

        net_charge.append(getNetCharge(x))
        Boman_index.append(getbomanindex(x))
        Aliphatic_Index.append(getAliphaticIndex(x))

        # 计算分子量
        molecular_weight = prot_param.molecular_weight()
        molecular.append(molecular_weight)

        # 计算等电点
        isoelectric_point = prot_param.isoelectric_point()
        charge_at_pH = prot_param.charge_at_pH(7.0)
        isoelectric.append(isoelectric_point)
        charge_pH.append(charge_at_pH)

        #芳香度
        aromaticity = prot_param.aromaticity()
        Aromaticity.append(aromaticity)

        # 计算疏水性值（gravy）
        gravy_value = prot_param.gravy()
        gravy.append(gravy_value)

        instability_index = prot_param.instability_index()
        instability.append(instability_index)

        secondary_structure_fraction = prot_param.secondary_structure_fraction()
        Helix.append(secondary_structure_fraction[0])
        Turn.append(secondary_structure_fraction[1])
        Sheet.append(secondary_structure_fraction[2])
        sequenceLength.append(len(x))

        flexibility = prot_param.flexibility()

        """hydrophobicity_values = {"A": 0.5, "C": 0.2, "D": -0.3, "E": -0.4, "F": 1.0, "G": 0.3, "H": -0.2, "I": 0.7, "K": -0.5, "L": 0.6, "M": 0.4, "N": -0.5, "P": 0.1, "Q": -0.5, "R": -0.7, "S": -0.2, "T": -0.1, "V": 0.6, "W": 0.9, "Y": 0.5}
        hydrophilicity_values = {"A": -0.5, "C": -0.2, "D": 0.3, "E": 0.4, "F": -1.0, "G": -0.3, "H": 0.2, "I": -0.7, "K": 0.5, "L": -0.6, "M": -0.4, "N": 0.5, "P": -0.1, "Q": 0.5, "R": 0.7, "S": 0.2, "T": 0.1, "V": -0.6, "W": -0.9, "Y": -0.5}
        """
        hydrophobicity_values = {
        'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 
        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 
        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 
        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }

        hydrophilicity_values = {
            'A': -0.5, 'R': 3, 'N': 0.2, 'D': 3, 'C': -1, 
            'Q': 0.2, 'E': 3, 'G': 0, 'H': -0.5, 'I': -1.8, 
            'L': -1.8, 'K': 3, 'M': -1.3, 'F': -2.5, 'P': 0, 
            'S': 0.3, 'T': -0.4, 'W': -3.4, 'Y': -2.3, 'V': -1.5
        }

        amphipathicity_index = calculate_amphipathicity(x, hydrophobicity_values, hydrophilicity_values)
        amphipathicity.append(amphipathicity_index)

        prot_param = peptides.Peptide(x)

        # 計算序列的自協方差指數
        auto_covariance = prot_param.auto_covariance(table)
        covariance.append(auto_covariance)

        # 計算序列的自相關性指數
        auto_correlation = prot_param.auto_correlation(table)
        correlation.append(auto_correlation)

        # 計算序列的疏水性矩
        hydrophobic_moment_value = prot_param.hydrophobic_moment()
        hydrophobic_moment.append(hydrophobic_moment_value)

        # 計算序列的質量值
        mass_value = prot_param.mass_shift()
        mass.append(mass_value)

        # 計算序列的膜位置分布
        membrane_position_profile_value = prot_param.membrane_position_profile()
        membrane_position_profile.append(membrane_position_profile_value)

        # 計算序列的質譜值
        mz_value = prot_param.mz()
        mz.append(mz_value)

        print("序列:", x)
        print("Amphipathicity Index:", amphipathicity_index)
        print("彈性:", flexibility)
        print("氨基酸组成：", aa_composition)
        print("分子量：", molecular_weight)
        print("pI:", isoelectric_point)
        print("pH7下的電荷:", charge_at_pH)
        print("疏水性值:", gravy_value)
        print("芳香度：", aromaticity)
        print("不穏定:", instability_index)
        print("二級:", secondary_structure_fraction)
        print("淨電荷:",net_charge[-1])

    new_data = pd.DataFrame({
        'sequence': data['sequence'],
        'gravy': gravy,
        'instability_index': instability,
        'Aliphatic_Index': Aliphatic_Index,
        'sequenceLength': sequenceLength,
        'isoelectric_point': isoelectric,
        'net_charge': net_charge,
        'molecular_weight': molecular,
        'charge_at_pH': charge_pH,
        'aromaticity': Aromaticity,
        'secondary_structure_fraction_Helix': Helix,
        'secondary_structure_fraction_Turn': Turn,
        'secondary_structure_fraction_Sheet': Sheet,
        'Boman_Index': Boman_index,
        'amphipathicity': amphipathicity,
        'correlation': correlation,
        'covariance': covariance,
        'hydrophobic_moment': hydrophobic_moment,
        'mass': mass,
        #'membrane_position_profile': membrane_position_profile,
        'mz': mz,
        
    })

    df = pd.DataFrame([ peptides.Peptide(s).descriptors() for s in data['sequence'] ])
    new_data = pd.concat([new_data, df], axis=1)
    new_data['MIC'] = data['MIC']
    new_data.to_csv(f'D:\\張\\科技部\\{bacteria}\\biopython-{bacteria}.csv', index=False, encoding='utf-8')