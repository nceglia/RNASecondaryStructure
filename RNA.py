"""
Classes:
Defining Secondary Structure Configuration
"""
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
    domain = []
    for i in xrange(len(sequence)):
        domain.append([])
        for j in xrange(len(sequence)):
            if i == j:
                domain[i].append([0])
            elif sequence[i]==sequence[j]:
                domain[i].append([0])
            else:
                domain[i].append([0,1])
    return domain

def CostStructure(structure):
    cost =0
    for i in xrange(len(structure)):
        count = 0
        for j in xrange(len(structure[i])):
            count += structure[i][j]
        if count < 2:
            cost += count
        else:
            cost += (1-count)
    return cost + stack(structure)

def TestCost(graph):
    cost = 0
    prev = graph[0]
    for i in graph[1:]:
        if (prev == i):
            cost = cost + 1
        prev = i
    return cost

