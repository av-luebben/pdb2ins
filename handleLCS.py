__author__= 'anna'
"""
Bridge between pdb2ins and longest common substring.
Maybe make some pretty data representation?
"""
# import time


class HandleLCS(object):
    def __init__(self):
        self.naturalAA = []
        self.naturalAAfrom3to1 = {}
        self.chainDict = {}
        self.shortAAChainDict = {}

    naturalAA = ['Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met',
                 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val']
    naturalAAfrom3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P',
                         'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R',
                         'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    def convert3to1(self):
        """
        take dict containing chains & 3 digit aa code and convert to string containing 1 digit aa code.
        3 digit dictionary: 'chainDict', Key = chainID, value = list of 3 digit AA,
        1 digit dictionary: 'shortAAChainDict', key = chainID, value = string of 1 digit AA, X as dummy for non natural
        AA, water residues are skipped.
        :return:
        """
        # dict = Header.getResiDict()
        for key in self.chainDict.keys():
            newAAList = []
            # print key, len(self.chainDict[key]), self.chainDict[key]
            for aa in self.chainDict[key]:
                if aa.capitalize() in self.naturalAA:
                    oneDigitAA = self.naturalAAfrom3to1[aa]
                    newAAList.append(oneDigitAA)
                elif aa.capitalize() == 'HOH':
                    pass
                else:
                    newAAList.append('X')
            self.shortAAChainDict[key] = ''.join(newAAList)
        # for key in self.shortAAChainDict.keys():
        #     print key, len(self.shortAAChainDict[key]), self.shortAAChainDict[key]

    def getAllInfo(self, chainDict):
        """
        dict = Header.getResiDict()
        Get a dictionary of strings which are later given to the longest common substring module. Dictionary: key= chainID,
        value =string of chain in 1-letter-abbreviation of natural amino acids in the chain.
        :param self:
        :param chainDict: dictionary; key=chain, value=string
        :return:
        """
        self.chainDict = chainDict
        self.convert3to1()


# def main():
#     global start_time
#     print 'INFO: +++ starting longest common substring search. +++ '
#     start_time = time.time()
#     main.__doc__ = head
#
#     print main.__doc__
#     #if version[:3] < '2.7':
#     #    print 'Please consider upgrading your python. This program requires python 2.7'
#     #    exit()
#     Data()
#
# if __name__ == "__main__":
#     main()
