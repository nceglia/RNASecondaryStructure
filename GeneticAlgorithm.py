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
import argparse

class GeneticAlgorithm(object):
    def __init__(self,domain,fitness,crossover=0.5,mutation=0.1,popSize=100,maxIter=600,maxUnchanged=500,tol=0.00001):
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
        self.count = 0
        self.prob = []
        self.scores = []
        self.contig = 0
        self.max_contig = maxUnchanged
        self.tol =tol

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
        if self.contig > self.max_contig:
            stopping = True
            print("No Change in {0} Iterations!".format(self.max_contig))
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

    def mutate(self,string,window=5):
        dim = numpy.array(string).shape
        mutatedString = numpy.empty(dim)
        for i in xrange(len(string)):
            for j in xrange(len(string[i])):
                if random.uniform(0.0, 1.0) <= self.mutation:
                    if len(self.domain[i][j]) == 2:
                        bindings = 0
                        window_size = 0
                        for l in xrange(-window,window,1):
                            for z in xrange(-window,window,1):
                                try:
                                    if string[i+l][j+z] == 1.0:
                                        bindings+=1
                                    window_size +=1
                                except Exception:
                                    continue
                        bound = float(bindings)/(window_size)
                        unbound = 1.0 - (bound)
                        choice = numpy.random.choice(self.domain[i][j],p=[unbound,bound])
                    else:
                        choice = 0.0
                    mutatedString[i][j] = choice
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
            bestChild = self.population[numpy.argmin(scores)]
            self.fittestChild = bestChild
        self.avgScore = numpy.mean(scores)
        print "Generation {2} Current Minimum: {0}, Average Score {1}".format(self.bestScore,self.avgScore,self.currentGeneration),
        return self.bestScore,self.minScore, self.avgScore, max(scores)

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
        self.run_times = []
        print("Initializing Population...")
        start_time = time.time()
        self.initialPopulation()
        self.evaluate()
        print "\t\tElapsed Time: ",time.time() - start_time, "Seconds"
        last_best = self.bestScore
        self.contig = 0
        self.score_watch = []
        while not self.stop():
            start_time = time.time()
            self.nextGeneration()
            self.score_watch.append(self.evaluate())
            if last_best - self.bestScore < self.tol:
                self.contig += 1
            else:
                self.contig = 0
            elapsed = time.time() - start_time
            print "\t\tElapsed Time: ",elapsed, "Seconds"
            self.run_times.append(elapsed)
            last_best = self.bestScore
        return self.bestScore, self.fittestChild

def runit(test,crossover=0.7,mutation=0.1,iterations=600,unchanged=50,draw=False):
    testcase = args.test
    print("**PREDICTING RNA SECONDARY STRUCTURE**")
    data = Data.Data("{0}.txt".format(testcase))
    seq = data.getSequence()
    algorithm = GeneticAlgorithm(StructureDomain(seq),
                                 lambda x: -1.0*CostStructure(x),
                                 crossover=args.crossover,
                                 mutation=args.mutation,
                                 maxIter=args.iterations,
                                 maxUnchanged=args.unchanged)
    bestScore, structure = algorithm.run()
    run_times = algorithm.run_times
    score_watch = algorith.score_watch
    print("Optimization Complete!")
    best_perc = data.testStructure(structure)
    if draw:
        data.drawDotBracket(structure,"{0}_{1}_Output".format(testcase.capitalize(),i))
    return run_times, score_watch, best_perc


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-t', dest='test', type=str, help='Test Case to Run')
    parser.add_argument('-d', dest='draw',help='Draw structure',action='store_true')
    parser.add_argument('-m',type=float,dest='mutation',default=0.1)
    parser.add_argument('-c',type=float,dest='crossover',default=0.7)
    parser.add_argument('-i',type=int,dest='iterations',default=600)
    parser.add_argument('-u',type=int,dest='unchanged',default=50)
    args = parser.parse_args()
    runit(args.test,args.crossover,args.mutation,args.iterations,args.unchanged,args.draw)
