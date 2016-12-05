"""
Load Fasta Inputs and Test Cases Samples
"""
import random
import numpy
import subprocess
import os
import time
from subprocess import Popen, PIPE, STDOUT
# seq1="ACUUCGCAAGCACGCGUAGGGAAAGGCACCAUGUAUCACGAUAUUACAUACUAAGAGCGU\
# CAACGUGAAUACCUGCUGGAUACUGUGUGGGCCGUGGUGAAAGUUUGAUCCGCAAAGCAG\
# CCCCUGUAACUGUACUCGCGGCAAGAGCAUCGCAGCAGUAUGUGCGUCUGAAUGCGACAC\
# GGAAGGCACGGCGGGACCCA"
# seq1 = numpy.array(list(seq1))
# testSq = [random.randint(0, 10) for i in range(100)]

class Data(object):
    def __init__(self,ifile):
        rows = open(ifile,"r").read().splitlines()
        self.sequence = numpy.array(list(rows[0]))
        self.solutions = []
        for row in rows[1:]:
            self.solutions.append(eval(row))

    def getSequence(self):
        return self.sequence

    @staticmethod
    def pairing(structure):
        result = []
        for i in xrange(len(structure)):
            for j in xrange(len(structure)):
                if float(structure[i][j]) == 1.0:
                    if (j,i) not in result:
                        result.append((i,j))
        for i,res in enumerate(result):
            if res[0] > res[1]:
                result[i] = (res[1],res[0])
        results = list(set(result))
        return results

    def testStructure(self,structure):
        results = self.pairing(structure)
        print "Returned Pairs: ", " ".join(map(str,results))
        for i, solution in enumerate(self.solutions):
            hits = 0
            for pair in results:
                if pair in solution:
                    hits+=1
            if hits == len(solution):
                print "{0}% Matched with Solution {1}".format(float(hits)/len(solution)*100.0,i)
        return

    def convertDotBracket(self,structure):
        dotbracket = ["." for _ in self.sequence]
        pairs = self.pairing(structure)
        for pair in pairs:
            dotbracket[pair[0]] = "("
            dotbracket[pair[1]] = ")"
        return "".join(dotbracket)

    def drawDotBracket(self,structure,prefix):
        pipe_data = ""
        pipe_data += ">{0}\n".format(prefix)
        pipe_data += "".join(list(self.sequence))+"\n"
        pipe_data += "".join(self.convertDotBracket(structure)) + "\n"
        p = Popen(['RNAPlot'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p.communicate(input=pipe_data)
