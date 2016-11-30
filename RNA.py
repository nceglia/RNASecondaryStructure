"""
Classes:
Defining Secondary Structure Configuration
Computing Energy of Secondary Structure
"""

def StructureDomain(sequence):
    """
    Will Add Structure Dependent Constraints
    0 = Stacked Pair
    1 = Bulge Loop
    2 = Interior Loop
    3 = K-Multiloop
    """
    #return [xrange(4) for _ in xrange(size)]
    return [range(10) for _ in range(100)]

def TestCost(graph):
    cost = 0
    prev = graph[0]
    for i in graph[1:]:
        if (prev == i):
            cost = cost + 1
        prev = i
    return cost

def GibbsFreeEnergy(configuration):
    totalEnergy = 0
    return totalEnergy
