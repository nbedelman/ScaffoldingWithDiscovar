from Bio.Seq import Seq
from Bio import SeqIO
import os
import itertools
import string
import copy
from operator import itemgetter
from ErrorClasses import *
from ContigClass import *
from ScaffoldClass import *

class Part(object):
    def __init__(self, partTuple):
        self.backbone=partTuple[0]

        self.type = type(self.backbone)
        self.name = self.backbone.getName()
        self.backboneLength = self.backbone.getLength()
        self.score = self.backbone.getScore()
        self.strand = partTuple[3]
        if self.type == Scaffold:
            self.originalStrand = self.backbone.getOriginalStrand()
        else:
            self.originalStrand = "NA"
        self.start = partTuple[1]
        self.end = partTuple[2]
    def getBackbone(self):
        return self.backbone
    def getType(self):
        if self.type == Segment:
            return Contig
        return self.type
    def getName(self):
        return self.name
    def getBackboneLength(self):
        return self.backboneLength
    def getScore(self):
        return self.score
    def getStrand(self):
        return self.strand
    def getOriginalStrand(self):
        return self.originalStrand
    def getStart(self):
        return self.start
    def getEnd(self):
        return self.end
    def getLength(self):
        return self.getEnd()- self.getStart()
    def printPart(self):
        output="%s,%i,%i,%s\n" % (self.getName(), self.getStart(), self.getEnd(), self.getStrand())
        return output
    def exportPart(self):
        return (self.getBackbone(), self.getStart(), self.getEnd(), self.getStrand())
    def flipPart(self):
        if self.getStrand() == '+':
            self.strand = '-'
        else:
            self.strand = '+'
        newEnd = self.getBackboneLength() - self.getStart()
        newStart = self.getBackboneLength() - self.getEnd()
        self.end = newEnd
        self.start = newStart
