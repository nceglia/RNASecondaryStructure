"""
Main Algorithm Class
"""
from RNA import StructureDomain,GibbsFreeEnergy
from Data import Data
import random
import logger
import numpy
import sys

class GeneticAlgorithm(object):
    def __init__(self,domain,fitness,crossover=0.7,mutation=0.01,popSize=1000,maxIter=100,avgChange=None):
        self.mutation = mutation
        self.crossover = crossover
        self.popSize = popSize
        self.fitness = fitness
        self.domain = domain
        self.maxIter = maxIter
        self.currentGeneration = 0
        self.bestScore = sys.maxfloat
        self.avgScore =sys.maxfloat
        self.oldScore = sys.maxfloat
        self.avgChange = avgChange
        count = 0
        
    def initialPopulation(self):
        self.population = [self.randomChromosome(self.domain) for _ in xrange(self.popSize)]

    def randomChromosome(self,domain):
        string = [random.choice(domain) for x in self.domain]
        return string

    def stop(self):
        stopping = False
        if self.avgScore - self.oldScore < self.avgChange:
            count += 1
        if count == 5:
            stopping = True
            self.oldScore = self.avgScore
        if self.currentGeneration > maxIter:
            stopping = True
            logging.info("Maximum Iterations Reached!")
        return stopping

    def scoreFitness(self):
        scores = []
        for chrom in self.population:
            scores.append(self.fitness(chrom))
        return scores

    def probability(self,scores)
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
                 mutatedString.append(random.choice(self.domain[i]))
            else:
                mutatedString.append(loc)
        return mutatedString

    def splice(self,parentA,parentB,points=1):
        splicedString = []
        for _ in xrange(points):
            splicePoint = random.randint(0, len(self.domain))
            for index in xrange(0,splicePoint):
                splicedString[index] = parentA[index]
            for index in xrange(splicePoint,len(self.domain)):
                splicedString[index] = parentB[index]
        return splicedString

    def evaluate(self):
        scores = self.scoreFitness()
        probability = self.probability(scores)
        minScore = min(scores)
        if minScore <= self.bestScore:
            self.bestScore = minScore
            self.fittestChild = self.population[numpy.argmin(scores)]
        self.avgScore = mean(scores)
        logger.info("Current Minimum: {0}, Average Score {1}".format(self.bestScore,self.avgScore))

    def nextGeneration(self):
        logger.info("Crossing Over and Mutated Generation {0}".format(self.currentGeneration))
        nextPop = []
        for _ in self.popSize:
            parentA, parentB = numpy.random.choice(nextPop, 2, p=prob)
            if random.uniform(0.0, 1.0) <= self.crossover:
                child = self.splice(parentA, parentB,points=1)
            else:
                random.choice([parentA,parentB])
            nextPop.append(self.mutate(child))
        self.currentGeneration += 1
        self.population = nextPop

    def run(self):
        logger.info("Initializing Population...")
        self.initialPopulation()
        self.evaluate()
        while not self.stop():
            self.nextGeneration()
            self.evaluate()
        return self.bestScore, self.fittestChild


def main():
    logger.info("**PREDICTING RNA SECONDARY STRUCTURE**")
    testSeq = Data.getSequence()
    testDomain = StructureDomain(testSeq)
    algorithm = GeneticAlgorithm(testDomain,GibbsFreeEnergy)
    algorithm.run()
    logger.info("Optimization Complete!")
    #print(testSeq.structure)

if __name__ == '__main__':
    main()
