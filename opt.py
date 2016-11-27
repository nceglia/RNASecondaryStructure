
# coding: utf-8

# In[1]:

import random
import numpy
import sys
import logging
testSq = [random.randint(0, 5) for i in range(50)]
class Data(object):
    def __init__(self,source):
        self.source = source

    def getSequence(self):
        return testSq
def StructureDomain(sequence):
    """
    Will Add Structure Dependent Constraints
    0 = Stacked Pair
    1 = Bulge Loop
    2 = Interior Loop
    3 = K-Multiloop
    """
    return [xrange(5) for _ in xrange(50)]

def TestCost(graph):
    cost = 0
    prev = graph[0]
    for i in graph[1:]:
        if (prev == i):
            cost = cost + 1
        prev = i
    return cost
   


# In[119]:


class GeneticAlgorithm(object):
    def __init__(self,domain,fitness,crossover=0.7,mutation=0.01,popSize=1000,maxIter=100,avgChange=0.01):
        self.mutation = mutation
        self.crossover = crossover
        self.popSize = popSize
        self.fitness = fitness
        self.domain = domain
        self.maxIter = maxIter
        self.currentGeneration = 0
        self.bestScore = sys.float_info.max
        self.avgScore =sys.float_info.max
        self.oldScore = sys.float_info.max
        self.avgChange = avgChange
        self.count = 0
        self.prob = []

    def initialPopulation(self):
        self.population = []
        for i in xrange(popSize):
            self.population.append(self.randomChromosome(self.domain))

    def randomChromosome(self,domain):
        string = [random.choice(self.domain[i]) for x in self.domain]
        return string

    def stop(self):
        stopping = False
        if self.avgScore - self.oldScore < self.avgChange:
            self.count += 1
        if self.count == 5:
            stopping = True
            self.oldScore = self.avgScore
        if self.currentGeneration > self.maxIter:
            stopping = True
            logging.info("Maximum Iterations Reached!")
        return stopping

    def scoreFitness(self):
        scores = []
        for chrom in self.population:
            scores.append(self.fitness(chrom))
        return scores

    def probability(self,scores):
        tmp = numpy.max(scores)
        scores -= tmp
        scores = numpy.exp(scores)
        tmp = numpy.sum(scores)
        scores /= tmp
        scores = 1.0 - scores
        return scores

    def mutate(self,string):
        mutatedString = []
        for loc in string:
            if random.uniform(0.0, 1.0) <= self.mutation:
                 mutatedString.append(random.choice(self.domain[0]))
            else:
                mutatedString.append(loc)
        return mutatedString

    def splice(self,parentA,parentB,points=1):
        splicedString = [None] * len(self.domain)
        for _ in xrange(points):
            splicePoint = random.randint(0, len(self.domain))
            for index in xrange(0,splicePoint):
                splicedString[index] = parentA[index]
            for index in xrange(splicePoint,len(self.domain)):
                splicedString[index] = parentB[index]
        return splicedString

    def evaluate(self):
        scores = self.scoreFitness()
        self.prob = self.probability(scores)
        minScore = min(scores)
        if minScore <= self.bestScore:
            self.bestScore = minScore
            self.fittestChild = self.population[numpy.argmin(scores)]
        self.avgScore = numpy.mean(scores)
        logging.info("Current Minimum: {0}, Average Score {1}".format(self.bestScore,self.avgScore))

    def nextGeneration(self):
        logging.info("Crossing Over and Mutated Generation {0}".format(self.currentGeneration))
        nextPop = []
        for _ in range(self.popSize):
            plength = len(self.population)
            index = list(range(plength))
            Aindex, Bindex = numpy.random.choice(index, 2, replace = False)
            parentA = self.population[Aindex]
            parentB = self.population[Bindex]
            if random.uniform(0.0, 1.0) <= self.crossover:
                child = self.splice(parentA, parentB,points=1)
            else:
                child = random.choice([parentA,parentB])
            nextPop.append(self.mutate(child))
        self.currentGeneration += 1
        self.population = nextPop

    def run(self):
        logging.info("Initializing Population...")
        self.initialPopulation()
        self.evaluate()
        while not self.stop():
            self.nextGeneration()
            self.evaluate()
        return self.bestScore, self.fittestChild


def main():
    logging.info("**PREDICTING RNA SECONDARY STRUCTURE**")
    testSeq = Data(testSq).getSequence()
    testDomain = StructureDomain(testSeq)
    algorithm = GeneticAlgorithm(testDomain,TestCost)
    algorithm.run()
    logging.info("Optimization Complete!")
    #print(testSeq.structure)




# In[121]:

a = main()


# In[ ]:



