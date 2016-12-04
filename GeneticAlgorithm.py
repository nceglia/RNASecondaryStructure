"""
Main Algorithm Class
"""
from RNA import StructureDomain,CostStructure
import Data
import random
import logging
import numpy
import sys
import math
import time

class GeneticAlgorithm(object):
    def __init__(self,domain,fitness,crossover=0.7,mutation=0.1,popSize=100,maxIter=500,avgChange=None):
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
        if self.currentGeneration > self.maxIter:
           stopping = True
           print("Maximum Iterations Reached!")
        return stopping

    def scoreFitness(self):
        scores = []
        for i, chrom in enumerate(self.population):
            scores.append(self.fitness(chrom))
        self.scores = scores
        return scores

    def probability(self,scores):
        norm = [1.0/(max(scores)-min(scores))*(float(x)-max(scores))+1.0 for x in scores]
        e = numpy.exp(-numpy.array(norm) / 1.0)
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

    def splice(self,parentA,parentB):
        dimBefore = numpy.array(parentA).shape
        flattenedA = numpy.array(parentA).flatten()
        flattenedB = numpy.array(parentB).flatten()
        splicedString = [None] * len(flattenedA)
        splicePoint = random.randint(0, len(flattenedA))
        for index in xrange(0,splicePoint):
            splicedString[index] = flattenedA[index]
        for index in xrange(splicePoint,len(flattenedA)):
            splicedString[index] = flattenedB[index]
        splicedString = numpy.array(splicedString).reshape(dimBefore)
        return splicedString

    def evaluate(self):
        scores = self.scoreFitness()
        self.prob = self.probability(scores)
        minScore = min(scores)
        if minScore <= self.bestScore:
            self.bestScore = minScore
            self.fittestChild = self.population[numpy.argmin(scores)]
        self.avgScore = numpy.mean(scores)
        print "Current Minimum: {0}, Average Score {1}".format(self.bestScore,self.avgScore),

    def nextGeneration(self):
        nextPop = []
        for _ in xrange(self.popSize):
            plength = len(self.population)
            index = list(range(plength))
            Aindex, Bindex = numpy.random.choice(index, 2, replace = False, p=self.prob)
            parentA = self.population[Aindex]
            parentB = self.population[Bindex]
            if random.uniform(0.0, 1.0) <= self.crossover:
                child = self.splice(parentA, parentB)
            else:
                child = random.choice([parentA,parentB])
            nextPop.append(self.mutate(child))
        self.currentGeneration += 1
        self.population = nextPop

    def run(self):
        print("Initializing Population...")
        start_time = time.time()
        self.initialPopulation()
        self.evaluate()
        print "\t\tElapsed Time: ",time.time() - start_time, "Seconds"
        while not self.stop():
            start_time = time.time()
            self.nextGeneration()
            self.evaluate()
            print "\t\tElapsed Time: ",time.time() - start_time, "Seconds"
        return self.bestScore, self.fittestChild


def main():
    print("**PREDICTING RNA SECONDARY STRUCTURE**")
    data = Data.Data("test2.txt")
    seq = data.getSequence()
    #testDomain = StructureDomain(seq1)
    algorithm = GeneticAlgorithm(StructureDomain(seq),lambda x: -1.0*CostStructure(x))
    bestScore, structure = algorithm.run()
    print("Optimization Complete!")
    data.testStructure(structure)
if __name__ == '__main__':
    main()
