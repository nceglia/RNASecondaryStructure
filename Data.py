"""
Load Fasta Inputs and Test Cases Samples
"""
import random

seq1="ACUUCGCAAGCACGCGUAGGGAAAGGCACCAUGUAUCACGAUAUUACAUACUAAGAGCGU\
CAACGUGAAUACCUGCUGGAUACUGUGUGGGCCGUGGUGAAAGUUUGAUCCGCAAAGCAG\
CCCCUGUAACUGUACUCGCGGCAAGAGCAUCGCAGCAGUAUGUGCGUCUGAAUGCGACAC\
GGAAGGCACGGCGGGACCCA"

testSq = [random.randint(0, 10) for i in range(100)]


def Data(object):
    def __init__(self,source):
        self.source = source

    def getSequence(self):
        #return seq1
        return testSq



