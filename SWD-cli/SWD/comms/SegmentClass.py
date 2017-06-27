from Bio.Seq import Seq
from Bio import SeqIO
import os
import itertools
import string
import copy
from operator import itemgetter
from ScaffoldClass import *
from ContigClass import *
from ErrorClasses import *

class Segment(object):
    def __init__(self,segment):
        atts=segment.split("\t")
        self.chrom=atts[0]
        self.relStart=int(atts[1])
        self.relEnd=int(atts[2])
        self.name=atts[3]
        self.score=atts[4]
        self.strand=atts[5]
        self.start=int(atts[6])
        self.end=int(atts[7])
        self.color=atts[8].strip("\n")
        self.overlap=[]
        self.contig=''
        self.isFlipped=False
        self.originalStrand=self.strand
    def flipSegment(self):
        self.isFlipped=not self.isFlipped
        newOverlaps=[]
        for overlap in self.getOverlap():
            newOverlap=copy.copy(overlap)
            newOverlap.flipScaffold()
            newOverlaps.append(newOverlap)
        self.overlap=newOverlaps
    def getOriginalStrand(self):
        return self.originalStrand
    def getConLength(self):
        return self.getRelEnd()-self.getRelStart()
    def getContig(self):
        return self.contig
    def getConStart(self):
        if self.getOriginalStrand() == '-':
            return self.getRelEnd()-self.getEnd()
        else:
            return self.getStart()-self.getRelStart()
    def getConEnd(self):
        if self.getOriginalStrand() == '-':
            return self.getConLength()-(self.getStart()-self.getRelStart())
        else:
            return self.getConLength()-(self.getRelEnd()-self.getEnd())
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
        if self.isFlipped:
            if self.strand == '+':
                return '-'
            else:
                return '+'
        else:
            return self.strand
    def getRelStart(self):
        return self.relStart
    def getRelEnd(self):
        return self.relEnd
    def getScore(self):
        return self.score
    def getColor(self):
        return self.color
    def getScafStart(self):
        if len(self.getOverlap()) > 1:
            raise (str(self.getName())+ "overlaps at least two scaffolds...")
        else:
            for scaffold in self.getOverlap():
                if not scaffold.isFlipped:
                    return self.getStart()-scaffold.getStart()
                else:
                    return scaffold.getEnd()-self.getEnd()
    def getScafEnd(self):
        if len(self.getOverlap()) > 1:
            raise str(self.getName())+ "overlaps at least two scaffolds..."
        else:
            for scaffold in self.getOverlap():
                if not scaffold.isFlipped:
                    return scaffold.getLength()-(scaffold.getEnd()-self.getEnd())
                else:
                    return scaffold.getLength()-(self.getStart()-scaffold.getStart()) +1
    def getConStartPos(self):
        if not self.isFlipped:
            return self.getStart()-self.getRelStart()
        else:
            return self.getRelEnd()-self.getEnd() +1
    def getConEndPos(self):
        if not self.isFlipped:
            return self.getConLength() - (self.getRelEnd()-self.getEnd())
        else:
            return self.getConLength() - (self.getStart()-self.getRelStart()) +1

    def getDistanceFromScafStart(self):
        return self.getScafStart()

    def getDistanceFromScafEnd(self):
        for i in self.getOverlap():
            return i.getLength()-self.getScafEnd()

    def findOverlaps(self, scaffoldList, markScaffolds):
        overlaps=[]
        for scaffold in scaffoldList:
            if scaffold.getChrom() == self.getChrom():
                oLength=0
                if (self.getStart() >= scaffold.getStart() and self.getStart() <= scaffold.getEnd()):
                    if (self.getEnd() >= scaffold.getStart() and self.getEnd() <= scaffold.getEnd()):
                        oLength=self.getLength()
                        overlaps.append((scaffold, oLength))
                        if (not self.checkIfUsed(scaffold.getOverlaps(), self.getContig())) and markScaffolds:
                            scaffold.overlaps.append(self.getContig())
                    else:
                        oLength=scaffold.getEnd()-self.getStart()
                        if (not self.checkIfUsed(scaffold.getOverlaps(), self.getContig())) and markScaffolds:
                            scaffold.overlaps.append(self.getContig())
                elif (self.getEnd() >= scaffold.getStart() and self.getEnd() <= scaffold.getEnd()):
                    oLength=self.getEnd()-scaffold.getStart()
                    overlaps.append((scaffold, oLength))
                    if (not self.checkIfUsed(scaffold.getOverlaps(), self.getContig())) and markScaffolds:
                        scaffold.overlaps.append(self.getContig())
        if len(overlaps) > 1:
            ordered=sorted(overlaps, key=itemgetter(1), reverse=True)
            best=ordered[0]
            onlyOverlap=best[0]
            self.overlap=[onlyOverlap,]
            return [onlyOverlap,]
        else:
            overlap=overlaps[0]
            onlyOverlap=overlap[0]
            self.overlap=[onlyOverlap,]
            return [onlyOverlap,]

    def getOverlap(self):
        return self.overlap

    def checkIfUsed(self,aList, anObject):
        used=False
        for item in aList:
            if item.getName() == anObject.getName():
                used=True
                break
        return used
