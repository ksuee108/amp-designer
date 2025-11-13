from Bio.SeqUtils.ProtParam import ProteinAnalysis
import peptides
class Bio_analysis():
    """
    这个类用于计算蛋白质序列的各种特征和性质。
    """

    def __init__(self, seq):
        """
        初始化類並檢查序列是否為有效的蛋白質序列。

        :param sequence: 輸入的蛋白質序列 (字串)
        """
        self.seq = seq
        self.prot_param = ProteinAnalysis(seq)
        self.prot_param2 = peptides.Peptide(seq)
        # Boman index data
        self.boman = {
            "A": 1.81, "C": 1.28, "D": -8.72, "E": -6.81, "F": 2.98,
            "G": 0.94, "H": -4.66, "I": 4.92, "K": -5.55, "L": 4.92,
            "M": 2.35, "N": -6.64, "P": 0, "Q": -5.54, "R": -14.92,
            "S": -3.4, "T": -2.57, "V": 4.04, "W": 2.33, "Y": -0.14
        }
        self.table = peptides.tables.HYDROPHOBICITY["KyteDoolittle"]
 
    def get_flexibility(self):
        """
        獲取蛋白質的柔性指數 (Flexibility)。
        這是基於 B-values (實驗 B 因子) 的計算結果。

        :return: 柔性指數的數值列表，每個殘基對應一個值
        """
        return self.prot_param.flexibility()
 
    def get_aa_composition(self):
        """
        獲取氨基酸的百分比組成。

        :return: 包含氨基酸百分比組成的字典
        """
        return self.prot_param.get_amino_acids_percent()

    def get_molecular_weight(self):
        """
        計算蛋白質的分子量。

        :return: 蛋白質的分子量 (浮點數)
        """
        return self.prot_param.molecular_weight()

    def get_isoelectric_point(self):
        """
        獲取蛋白質的等電點。

        :return: 蛋白質的等電點 (浮點數)
        """
        return self.prot_param.isoelectric_point()

    def get_charge_at_pH(self, pH=7.0):
        """
        獲取指定 pH 值下的蛋白質電荷。

        :param pH: pH 值，默認為 7.0
        :return: 蛋白質在指定 pH 值下的電荷 (浮點數)
        """
        return self.prot_param.charge_at_pH(pH)

    def get_aromaticity(self):
        """
        獲取蛋白質的芳香度。

        :return: 蛋白質的芳香度 (浮點數)
        """
        return self.prot_param.aromaticity()

    def get_gravy(self):
        """
        獲取蛋白質的疏水性 (GRAVY)。
        >0 = 疏水性
        <0 = 親水性

        :return: 蛋白質的疏水性 (浮點數)
        """
        return self.prot_param.gravy()

    def get_instability_index(self):
        """
        獲取蛋白質的不穩定指數。

        :return: 不穩定指數 (浮點數)
        """
        return self.prot_param.instability_index()

    def get_secondary_structure_fraction(self):
        """
        獲取蛋白質的次級結構分數。

        :return: 包含 alpha 螺旋、beta 板和 beta 轉角百分比的元組
        """
        return self.prot_param.secondary_structure_fraction()

    def get_net_charge(self):
        """
        獲取蛋白質在 pH 7.0 下的淨電荷。

        :return: 蛋白質的淨電荷 (浮點數)
        """
        net_charge = 0.0
        for aa in self.seq:
            if aa in ["R", "K", "H"]:
                net_charge += 1.0
            elif aa in ["D", "E"]:
                net_charge -= 1.0
        return net_charge

    def get_aliphatic_index(self):
        """
        獲取蛋白質的脂肪族指數。

        :return: 脂肪族指數 (浮點數)
        """
        a_factor = 2.9
        b_factor = 3.9
        a_count = self.seq.count("A")
        v_count = self.seq.count("V")
        il_count = self.seq.count("I") + self.seq.count("L")
        
        total_length = len(self.seq)
        a_percent = (a_count / total_length) * 100
        v_percent = (v_count / total_length) * 100
        il_percent = ((il_count) / total_length) * 100
        
        return a_percent + (a_factor * v_percent) + (b_factor * il_percent)

    def get_boman_index(self):
        """
        獲取蛋白質的 Boman 指數。

        :return: Boman 指數 (浮點數)
        """
        boman_sum = 0
        for aa in self.seq:
            boman_sum += self.boman.get(aa, 0)  # 防止不存在的鍵
        
        return -boman_sum / len(self.seq)
    
    def get_sequenceLength(self):
        """
        獲取蛋白質序列的長度。

        :return: 序列的長度
        """
        return len(self.seq)
    
    def get_molar_extinction_coefficient(self):
        """
        獲取蛋白質序列的摩爾消光係數。

        :return: 摩爾消光係數
        """
        return self.prot_param.molar_extinction_coefficient()
    
    def get_auto_correlation(self):
        """
        獲取蛋白質序列的自動協方差。
        """
        return self.prot_param2.auto_correlation(self.table)
    
    def get_auto_covariance(self):
        """
        獲取蛋白質序列的自協方差。
        """
        return self.prot_param2.auto_covariance(self.table)
    
    def get_hydrophobic_moenet(self):
        """
        獲取蛋白質序列的疏水矩。
        """
        return self.prot_param2.hydrophobic_moment()
    
    def get_mass(self):
        """
        獲取蛋白質序列的質量。
        """
        return self.prot_param2.mass_shift()
    
    def get_mz(self):
        """
        獲取蛋白質序列的質譜值。
        """
        return self.prot_param2.mz()
    
    def get_distance_matrix(self):
        """
        獲取蛋白質序列的距離矩陣。
        """
        return self.prot_param2.descriptors()