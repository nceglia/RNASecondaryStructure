"""
Load Sequences and Test Cases Samples
"""
import random
import numpy
import subprocess
import os
import time
from subprocess import Popen, PIPE, STDOUT
import difflib

class Data(object):
    def __init__(self,ifile,testType="pairs"):
        rows = open(ifile,"r").read().splitlines()
        self.sequence = numpy.array(list(rows[0]))
        self.solutions = []
        self.testType = testType
        for row in rows[1:]:
            try:
                self.solutions.append(eval(row))
            except SyntaxError:
                self.solutions.append(row)
                self.testType = "dotbracket"

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
        if self.testType == "pairs":
            results = self.pairing(structure)
            best_scores = []
            print "Returned Pairs: ", " ".join(map(str,results))
            for i, solution in enumerate(self.solutions):
                hits = 0
                for pair in results:
                    if pair in solution:
                        hits+=1
                best_scores.append(float(hits)/len(solution))
                print "{0}% Matched with Solution {1}".format(float(hits)/len(solution)*100.0,i)
            best_score = max(best_scores)
        elif self.testType == "dotbracket":
            dotbracket = self.convertDotBracket(structure)
            print "Returned DotBracket {0}".format(dotbracket)
            for i, solution in enumerate(self.solutions):
                perc = difflib.SequenceMatcher(None,dotbracket,str(solution))
                print "{0}% Matched with Solution {1}".format(float(perc.ratio())*100.0,i)
            best_score = perc.ratio()
        return best_score

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
        p = Popen(['RNAPlot','-t','0'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
        p.communicate(input=pipe_data)
