"""
Main Algorithm Class
"""
from RNA import StructureDomain,CostStructure
from Data import *
import random
import logging
import numpy
import sys

class GeneticAlgorithm(object):
    def __init__(self,domain,fitness,crossover=0.7,mutation=0.01,popSize=25,maxIter=100,avgChange=None):
        self.mutation = mutation
        self.crossover = crossover
        self.popSize = popSize
        self.fitness = fitness
        self.domain = domain
        self.maxIter = maxIter
        self.currentGeneration = 0
        self.bestScore = sys.float_info.max
        self.avgScore = sys.float_info.max
        self.oldScore = sys.float_info.max
        self.avgChange = avgChange
        self.count = 0
        self.prob = []
        self.scores = []

    def initialPopulation(self):
        self.population = []
        for i in range(self.popSize):
            self.population.append(self.randomChromosome(self.domain))
    def randomChromosome(self,domain):
        string = []
        for line in self.domain:
            string.append([random.choice(x) for x in line])
        return string

    def stop(self):
        stopping = False

        #if self.avgScore - self.oldScore < self.avgChange:
         #   self.count += 1
        #if self.count == 5:
         #   stopping = True
          #  self.oldScore = self.avgScore
        if self.currentGeneration > self.maxIter:
           stopping = True
           print("Maximum Iterations Reached!")
        # if (self.bestScore == 0):
        #     stopping = True
        return stopping

    def scoreFitness(self):
        scores = []
        print "In Score Fitness.."
        for i, chrom in enumerate(self.population):
            scores.append(self.fitness(chrom))
        self.scores = scores
        return scores

    def probability(self,scores):
        npa = numpy.array
        e = numpy.exp(-npa(scores) / 1.0)
        dist = e / numpy.sum(e)
        return dist

    def mutate(self,string):
        dim = numpy.array(string).shape
        mutatedString = numpy.empty(dim)
        for i in xrange(len(string)):
            for j in xrange(len(string[i])):
                if random.uniform(0.0, 1.0) <= self.mutation:
                 mutatedString[i][j] = random.choice(self.domain[i][j])
                else:
                    mutatedString[i][j] = string[i][j]
        return mutatedString

    def splice(self,parentA,parentB,points=1):
        dimBefore = numpy.array(parentA).shape
        flattenedA = numpy.array(parentA).flatten()
        flattenedB = numpy.array(parentB).flatten()
        splicedString = [None] * len(flattenedA)
        for _ in xrange(points):
            splicePoint = random.randint(0, len(flattenedA))
            for index in xrange(0,splicePoint):
                splicedString[index] = flattenedA[index]
            for index in xrange(splicePoint,len(self.domain)):
                splicedString[index] = flattenedB[index]
        return numpy.array(splicedString).reshape(dimBefore)

    def evaluate(self):
        print "Scoring Fitness..."
        scores = self.scoreFitness()
        print "Scoring Probability"
        self.prob = self.probability(scores)
        minScore = min(scores)
        if minScore <= self.bestScore:
            self.bestScore = minScore
            self.fittestChild = self.population[numpy.argmin(scores)]
        self.avgScore = numpy.mean(scores)
        print("Current Minimum: {0}, Average Score {1}".format(self.bestScore,self.avgScore))

    def nextGeneration(self):
        print("Crossing Over and Mutated Generation {0}".format(self.currentGeneration))
        nextPop = []
        for _ in xrange(self.popSize):
            plength = len(self.population)
            index = list(range(plength))
            Aindex, Bindex = numpy.random.choice(index, 2, replace = False, p=self.prob)
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
        print("Initializing Population...")
        self.initialPopulation()
        print("Evaluating Initial Population..")
        self.evaluate()
        while not self.stop():
            self.nextGeneration()
            self.evaluate()
        return self.bestScore, self.fittestChild


def main():
    print("**PREDICTING RNA SECONDARY STRUCTURE**")
    # testSeq = Data.Data(seq1).getSequence()
    testDomain = StructureDomain(seq1)
    algorithm = GeneticAlgorithm(StructureDomain(seq1),lambda x: -1.0*CostStructure(x))
    algorithm.run()
    print("Optimization Complete!")

if __name__ == '__main__':
    main()
