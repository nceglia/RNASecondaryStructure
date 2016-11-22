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
    return [xrange(4) for _ in xrange(size)]


def GibbsFreeEnergy(configuration):
    totalEnergy = 0
    return totalEnergy