"""
Defining Secondary Structure Configuration
"""
import numpy
from Data import Data

def StructureDomain(sequence):
    domain = numpy.empty((len(sequence),len(sequence),2))
    for i in xrange(len(sequence)):
        for j in xrange(len(sequence)):
            if i == j:
                domain[i][j][0] = 0
                domain[i][j][1] = 0
            elif (sequence[i]=='A' and sequence[j]=='U') or (sequence[i]=='U' and sequence[j]=='A'):
                domain[i][j][0] = 0
                domain[i][j][1] = 1
            elif (sequence[i]=='G' and sequence[j]=='C') or (sequence[i]=='C' and sequence[j]=='G'):
                domain[i][j][0] = 0
                domain[i][j][1] = 1
            else:
                domain[i][j][0] = 0
                domain[i][j][1] = 0
    return domain

def Penalty(structure):
    pairs = Data.pairing(structure)
    bases = list(numpy.array(pairs).flatten())
    cost = 0
    for base in list(set(bases)):
        count = bases.count(base)
        if count > 1:
            cost += count
    return cost

def CostStructure(structure,weight=-10.0):
    basepairs = 0
    for row in structure:
        for item in row:
            if item==1.0:
                basepairs+=1
    penalty = weight *Penalty(structure)
    return basepairs + penalty

