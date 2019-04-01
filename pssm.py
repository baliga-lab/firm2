# A class designed to hold a position specific scoring matrix
# and be able to output this in many different formats
#
# Variables:
# name - name
# eValue - the significance of the motif
# sites - genes that have the motif
# genes - the genes that have the motif
# matrix - a matrix of x by 4 [A, C, G, T]
class pssm:
    # Initialize the pssm
    def __init__(self, name=None, sites=None, evalue=None, pssm=None, genes=None):
        self.name = name
        self.sites = sites
        self.evalue = evalue
        self.matrix = pssm
        self.genes = genes

    def consensus_motif(self, lim1=0.6, lim2=0.8, three=0):
        """
        Returns the consensus word for a motif
        """
        consensus = ''
        for i in range(len(self.matrix)):
            consensus += self.col_consensus(self.matrix, i, lim1, lim2, three)
        return consensus

    def col_consensus(self, pssm, i, lim1, lim2, three):
        """
        Returns the consensus letter for a motif column
        """
        two_base_l = ['Y','R','W','S','K','M']
        three_base_l = ['V','H','D','B']
        conLet = 'N'

        if float(pssm[i][0]) >= lim1:
            conLet = 'A'
        elif float(pssm[i][1]) >= lim1:
            conLet = 'C'
        elif float(pssm[i][2]) >= lim1:
            conLet = 'G'
        elif float(pssm[i][3]) >= lim1:
            conLet = 'T'
        else:
            two_base_c = [float(pssm[i][1]) + float(pssm[i][3]),
                          float(pssm[i][0]) + float(pssm[i][2]),
                          float(pssm[i][0]) + float(pssm[i][3]),
                          float(pssm[i][1]) + float(pssm[i][2]),
                          float(pssm[i][2]) + float(pssm[i][3]),
                          float(pssm[i][0]) + float(pssm[i][1])]
            three_base_c = [float(pssm[i][0]) + float(pssm[i][1]) + float(pssm[i][2]),
                            float(pssm[i][0]) + float(pssm[i][1]) + float(pssm[i][3]),
                            float(pssm[i][0]) + float(pssm[i][2]) + float(pssm[i][3]),
                            float(pssm[i][1]) + float(pssm[i][2]) + float(pssm[i][3])]
            pMax = 0
            for k in range(0,6):
                if two_base_c[k] > pMax:
                    pMax = two_base_c[k]
                    conLet = two_base_l[k]
            if not pMax > lim2 and three==1:
                for k in range(0,4):
                    if three_base_c[k] > pMax:
                        pMax = three_base_c[k]
                        conLet = three_base_l[k]
            if not pMax > lim2:
                conLet = 'N'
        return conLet
