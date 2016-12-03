"""
Classes:
Defining Secondary Structure Configuration
"""
import numpy
def TestStructureDomain(sequence):
    return [range(10) for _ in range(100)]

def stack(structure,weightFunc = lambda x: x-1):
    stackWeight = 0
    for i in xrange(len(structure)):
        for j in reversed(xrange(len(structure))):
            n = 0
            z = i
            p = j
            if structure[i][j] == 1:
                n = 1
            while True:
                if z+1 >= len(structure) or p-1 < 0:
                    break
                elif structure[z+1][p-1] != 1:
                    break
                else:
                    n += 1
                z += 1
                p -= 1
            if n > 0:
                stackWeight += weightFunc(n)
    return stackWeight

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

def CostStructure(structure,weight=6.0):
    cost =0
    # How many base pairs
    basepairs = 0
    for row in structure:
        for item in row:
            if item==1.0:
                basepairs+=1
    cost = basepairs
    penalty = 0
    for i in xrange(len(structure)):
        count = 0
        for j in xrange(len(structure[i])):
            if structure[i][j] == 1.0:
                count += 1
        if count > 1:
            penalty += weight * count
    return cost - penalty + stack(structure)

def TestCost(graph):
    cost = 0
    prev = graph[0]
    for i in graph[1:]:
        if (prev == i):
            cost = cost + 1
        prev = i
    return cost

