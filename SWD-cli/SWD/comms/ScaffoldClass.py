from Bio.Seq import Seq
from Bio import SeqIO
import os
import itertools
import string
import copy
from operator import itemgetter
from ContigClass import *
from SegmentClass import *
from ErrorClasses import *

class Scaffold(object):
    def __init__(self,entry):
        atts=entry.split("\t")
        self.chrom=atts[0]
        self.start=int(atts[1])
        self.end=int(atts[2])
        self.name=atts[3]
        self.score=int(atts[4])
        self.strand=atts[5].strip("\n")
        self.overlaps=[]
        self.isFlipped=False
    def getOriginalStrand(self):
        return self.strand
    def flipScaffold(self):
        self.isFlipped= not self.isFlipped
    def getOverlaps(self):
        return self.overlaps
    def getChrom(self):
        return self.chrom
    def getStart(self):
        return self.start
    def getEnd(self):
        return self.end
    def getName(self):
        return self.name
    def getLength(self):
        return self.getEnd()-self.getStart()
    def getStrand(self):
        if not self.isFlipped:
            return self.strand
        elif self.strand =='+':
            return '-'
        else:
            return '+'
    def getScore(self):
        return self.score
