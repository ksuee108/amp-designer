from Bio.SeqUtils.ProtParam import ProteinAnalysis
import peptides

class Bio_analysis():
    """
    This class calculates various properties and features of protein sequences.
    """

    def __init__(self, seq):
        """
        Initialize the class and validate the protein sequence.

        :param seq: Input protein sequence (string)
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
        Get the flexibility index of the protein sequence based on experimental B-values.

        :return: List of flexibility values for each residue
        """
        return self.prot_param.flexibility()

    def get_aa_composition(self):
        """
        Get the percentage composition of amino acids.

        :return: Dictionary containing amino acid percentages
        """
        return self.prot_param.get_amino_acids_percent()

    def get_molecular_weight(self):
        """
        Calculate the molecular weight of the protein.

        :return: Molecular weight (float)
        """
        return self.prot_param.molecular_weight()

    def get_isoelectric_point(self):
        """
        Get the isoelectric point (pI) of the protein.

        :return: Isoelectric point (float)
        """
        return self.prot_param.isoelectric_point()

    def get_charge_at_pH(self, pH=7.0):
        """
        Get the net charge of the protein at a specified pH.

        :param pH: pH value (default 7.0)
        :return: Net charge at the specified pH (float)
        """
        return self.prot_param.charge_at_pH(pH)

    def get_aromaticity(self):
        """
        Get the aromaticity of the protein sequence.

        :return: Aromaticity (float)
        """
        return self.prot_param.aromaticity()

    def get_gravy(self):
        """
        Get the Grand Average of Hydropathy (GRAVY) of the protein.
        Positive = hydrophobic, negative = hydrophilic

        :return: GRAVY value (float)
        """
        return self.prot_param.gravy()

    def get_instability_index(self):
        """
        Get the instability index of the protein.

        :return: Instability index (float)
        """
        return self.prot_param.instability_index()

    def get_secondary_structure_fraction(self):
        """
        Get the secondary structure fractions.

        :return: Tuple with fractions of alpha-helix, beta-sheet, and beta-turn
        """
        return self.prot_param.secondary_structure_fraction()

    def get_net_charge(self):
        """
        Calculate the net charge at pH 7.0.

        :return: Net charge (float)
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
        Calculate the aliphatic index of the protein.

        :return: Aliphatic index (float)
        """
        a_factor = 2.9
        b_factor = 3.9
        a_count = self.seq.count("A")
        v_count = self.seq.count("V")
        il_count = self.seq.count("I") + self.seq.count("L")
        
        total_length = len(self.seq)
        a_percent = (a_count / total_length) * 100
        v_percent = (v_count / total_length) * 100
        il_percent = (il_count / total_length) * 100
        
        return a_percent + (a_factor * v_percent) + (b_factor * il_percent)

    def get_boman_index(self):
        """
        Calculate the Boman index of the protein.

        :return: Boman index (float)
        """
        boman_sum = 0
        for aa in self.seq:
            boman_sum += self.boman.get(aa, 0)  # Default 0 if amino acid not found
        return -boman_sum / len(self.seq)
    
    def get_sequenceLength(self):
        """
        Get the length of the protein sequence.

        :return: Sequence length (int)
        """
        return len(self.seq)
    
    def get_molar_extinction_coefficient(self):
        """
        Get the molar extinction coefficient.

        :return: Molar extinction coefficient (float)
        """
        return self.prot_param.molar_extinction_coefficient()
    
    def get_auto_correlation(self):
        """
        Calculate the auto-correlation descriptors using the hydrophobicity table.
        """
        return self.prot_param2.auto_correlation(self.table)
    
    def get_auto_covariance(self):
        """
        Calculate the auto-covariance descriptors using the hydrophobicity table.
        """
        return self.prot_param2.auto_covariance(self.table)
    
    def get_hydrophobic_moenet(self):
        """
        Calculate the hydrophobic moment of the protein sequence.
        """
        return self.prot_param2.hydrophobic_moment()
    
    def get_mass(self):
        """
        Calculate the mass shift of the protein sequence.
        """
        return self.prot_param2.mass_shift()
    
    def get_mz(self):
        """
        Calculate the m/z values of the protein.
        """
        return self.prot_param2.mz()
    
    def get_distance_matrix(self):
        """
        Get the distance matrix of the protein sequence.
        """
        return self.prot_param2.descriptors()
    
    def get_amphipathicity(self):
        """
        Calculate the amphipathicity of the protein sequence.

        :return: Amphipathicity value (float)
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

        hydrophobicity_sum = sum(hydrophobicity_values[aa] for aa in self.seq)
        hydrophilicity_sum = sum(hydrophilicity_values[aa] for aa in self.seq)
        return (hydrophobicity_sum / len(self.seq)) - (hydrophilicity_sum / len(self.seq))
