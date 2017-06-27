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
    
    def uniqueBisectScafSearch(self, scafList):
        #Find a scaffold in which the entirety of the segment maps
        totScafs=len(scafList)
        trialIndex=totScafs/2
        if (self.getStart() >= scafList[trialIndex].getStart()) and (self.getEnd() <= scafList[trialIndex].getEnd()):
            return scafList[trialIndex]
        elif totScafs==1:
            raise NoUniqueScafError(self.getName())
        elif self.getStart() > scafList[trialIndex].getStart():
            return self.uniqueBisectScafSearch(scafList[trialIndex:])
        else:
            return self.uniqueBisectScafSearch(scafList[:trialIndex])
            
    def nonUniqueBisectScafSearch(self, scafList):
        #in cases where there is no unique match, find the scaffold or scaffolds that include at least part of the segment
        totScafs=len(scafList)
        trialIndex=totScafs/2
        outList=[]
        newSegs=[]
        #If the seg maps to the end of the scaffold, add it to the results, and see if it also maps to the next scaffold
        #also, cut off the segment at the scaffold edge.
        if (self.getStart() >= scafList[trialIndex].getStart()) and (self.getStart() <= scafList[trialIndex].getEnd()):
            startSeg=copy.copy(self)
            startSeg.end=scafList[trialIndex].getEnd()
            startSeg.name=self.getName()+"a"
            newSegs.append(startSeg)
            outList.append(scafList[trialIndex])
            try:
                if (self.getEnd() >= scafList[trialIndex+1].getStart()) and (self.getEnd() <= scafList[trialIndex+1].getEnd()):
                    endSeg=copy.copy(self)
                    endSeg.start=scafList[trialIndex+1].getStart()
                    endSeg.name=self.getName()+"b"
                    newSegs.append(endSeg)
            except IndexError:
                self.end=newSegs[0].end
                return outList
            if len(newSegs)==1:
                self.end=newSegs[0].end
                return outList
            else:
                return newSegs
        
        #If the seg maps to the beginning of the scaffold, add it to the results, and see if it also maps to the previous scaffold
        #also, cut off the segment at the scaffold edge.
        elif (self.getEnd() >= scafList[trialIndex].getStart()) and (self.getEnd() <= scafList[trialIndex].getEnd()):
            endSeg=copy.copy(self)
            endSeg.start=scafList[trialIndex].getStart()
            endSeg.name=self.getName()+"b"
            newSegs.append(endSeg)
            outList.append(scafList[trialIndex])
            try:
                if (self.getStart() >= scafList[trialIndex-1].getStart()) and (self.getStart() <= scafList[trialIndex-1].getEnd()):
                    startSeg=copy.copy(self)
                    startSeg.end=scafList[trialIndex-1].getEnd()
                    startSeg.name=self.getName()+"a"
                    newSegs.append(startSeg)
            except IndexError:
                self.start=newSegs[0].start
                return outList
            if len(newSegs)==1:
                self.start=newSegs[0].start
                return outList
            else:
                return newSegs
              
        elif totScafs==1:
            raise NoUniqueScafError(self.getName())
        elif self.getStart() > scafList[trialIndex].getStart():
            return self.nonUniqueBisectScafSearch(scafList[trialIndex-1:])
        else:
            return self.nonUniqueBisectScafSearch(scafList[:trialIndex+1])

    def findOverlaps(self, scaffoldDict, markScaffolds):
        #find the scaffold(s) to which this segment maps. 
        #updated 6/27/17 to use a dictionary structure for the scaffolds instead of a list.
        #if there is a single scaffold, return the scaffold. If there are multiple, return segments split at the 
        #scaffold break and run them through again
        overlaps=[]
        #get the non-N scaffolds that are on the same chromosome 
        sameChrom=scaffoldDict[self.getChrom()][0]
        #use bisection search to try to find a scaffold to which the entire segment aligns. This should be almost all of them.
        try:
            bestOverlap=self.uniqueBisectScafSearch(sameChrom)
            overlaps.append(bestOverlap)

        except NoUniqueScafError:
            nonUniqueMatch=self.nonUniqueBisectScafSearch(sameChrom)
            if len(nonUniqueMatch) ==1:
                overlaps+=nonUniqueMatch
                bestOverlap=overlaps[0]
            else:
                return nonUniqueMatch

        if markScaffolds:
            for overlap in overlaps:
                if not (self.checkIfUsed(overlap.getOverlaps(), self.getContig())):
                    overlap.overlaps.append(self.getContig())
                    
        self.overlap=[bestOverlap,]
        return [bestOverlap,]


    def getOverlap(self):
        return self.overlap

    def checkIfUsed(self,aList, anObject):
        used=False
        for item in aList:
            if item.getName() == anObject.getName():
                used=True
                break
        return used
