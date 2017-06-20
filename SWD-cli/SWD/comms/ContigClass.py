from Bio.Seq import Seq
from Bio import SeqIO
import os
import itertools
import string
import copy
from operator import itemgetter
from ScaffoldClass import *
from SegmentClass import *
from ErrorClasses import *

class Contig(object):
    ungroupedChrom=''
    def __init__(self,bedFile):
        f= open(bedFile,"r")
        atts=f.readline().split("\t")
        self.chrom=atts[0]
        self.relStart=int(atts[1])
        self.relEnd=int(atts[2])
        self.name=atts[3]
        self.score=atts[4]
        self.strand=atts[5]
        self.start=int(atts[6])
        self.end=int(atts[7])
        self.color=atts[8]
        self.allSegments=[]
        for otherLine in f:
            segment=Segment(otherLine)
            segment.contig=self
            self.allSegments.append(segment)
        self.goodSegments=[]
        self.combinedSegments=[]
        self.connectors=[]
        self.inOrEx="e"
    def getScore(self):
        return self.score
    def getChrom(self):
        return self.chrom
    def getStart(self):
        return self.start
    def getEnd(self):
        return self.end
    def getName(self):
        return self.name
    def getLength(self):
        return self.getRelEnd()-self.getRelStart()
    def getStrand(self):
        return self.strand
    def getRelStart(self):
        return self.relStart
    def getRelEnd(self):
        return self.relEnd
    def getColor(self):
        return self.color
    def getAllSegments(self):
        return self.allSegments
    def getGoodSegments(self):
        return self.goodSegments
    def getCombinedSegments(self):
        return self.combinedSegments
    def getInOrEx(self):
        return self.inOrEx
    def cullSegments(self):
        #Tried to make it so even if the first match was to the ungroupedChrom, we could use it. That proved to cause more problems than it solved, so took that functionality out.
        #5/31/17
        chrom=self.getChrom()
        for segment in self.getAllSegments():
            if segment.getChrom() == chrom or segment.getChrom() == self.ungroupedChrom:
                duplicateSeg=False
                for goodSeg in self.getGoodSegments():
                    if goodSeg.getName() == segment.getName():
                        duplicateSeg=True
                        break
                if (duplicateSeg == False):
                    self.goodSegments.append(segment)

    def getConnectors(self):
        return self.connectors

    def combineSegments(self):
        '''here, combine segments that are very close to each other and really represent a single long mapping.
        Also get rid of very small (and therefore possibly spurious) alignments. Here, "small" is <1000 bp AND < 75% of the scaffold to which it maps'''
        if self.getGoodSegments() == []:
            raise NoSegmentError("Remember to cull segments before combining!")
        #elif len(self.getGoodSegments()) ==1:
        #    raise NotInformativeError("After culling, this contig no longer maps to two scaffolds")
        else:
            orderedSegs=self.orderSegs(self.getGoodSegments())
            combinedSegs=self.combineLists(orderedSegs)
            #if len(combinedSegs) == 1:
            #    raise NotInformativeError("After combining, this contig no longer maps to two scaffolds")
            #else:
            longEnough=[]
            for seg in combinedSegs:
                overlapLength=seg.getOverlap()[0].getLength()
                mappingPercent=float(seg.getConEnd()-seg.getConStart())/overlapLength
                if (seg.getConEnd()-seg.getConStart() > 1000) or (mappingPercent > .75) :
                    longEnough.append(seg)
            self.combinedSegments=longEnough
    def outputBed(self, segList, fileName):
        with open(fileName, "w") as o:
            for segment in segList:
                outList=[]
                outString=''
                for i in [segment.getChrom(),segment.getRelStart(),segment.getRelEnd(), segment.getName(), segment.getScore()\
                ,segment.getStrand(), segment.getStart(), segment.getEnd()]:
                    outString+=str(i)+"\t"
                outString+=segment.getColor()+"\n"
                o.write(outString)
        o.close()
    def findConnectors(self, scaffoldList, segType):
        connects=[]
        if segType=='good':
            segList=self.getGoodSegments()
            markScaffolds=False
        elif segType=='combined':
            segList=self.getCombinedSegments()
            markScaffolds=True
        for segment in segList:
            overlapper=segment.findOverlaps(scaffoldList, markScaffolds)
            for scaffold in overlapper:
                if (not scaffold in connects) and (not scaffold.getName() == "Ns"):
                    connects.append(scaffold)
        self.connectors=connects

    def orderSegs(self,segmentList):
        unordered=[]
        ordered=[]
        for i in segmentList:
            unordered.append((i,i.getStart()))
        ordered=sorted(unordered, key=itemgetter(1))
        segsOnly=[entry[0] for entry in ordered]
        return segsOnly



    def combineLists(self, goodSegs):
        combined=[]
        goodIndices=[]
        for i in (range(len(goodSegs)-1)):
            j=goodSegs[i]
            k=goodSegs[i+1]
            out = self.tryCombining(j,k)
            if out:
                goodIndices+=i,i+1
                combined.append(out)
        for l in range(len(goodSegs)):
            if not l in goodIndices:
                combined.append(goodSegs[l])
        if len(goodSegs) == len(combined):
            return goodSegs
        else:
            return self.combineLists(combined)


    def tryCombining(self, segOne, segTwo):
        if segOne.getStrand()!=segTwo.getStrand():
            return False
        elif segOne.getOverlap()[0].getName() != segTwo.getOverlap()[0].getName():
            return False
        elif (abs(segOne.getEnd()-segTwo.getStart()) < 1000 and self.closeOnContig(segOne, segTwo)) or ((segOne.getEnd() > segTwo.getStart())and (segOne.getStart() < segTwo.getStart())):
            newSegment=copy.copy(segOne)
            newSegment.relEnd=segTwo.getRelEnd()
            newSegment.end=segTwo.getEnd()
            return newSegment
        elif (abs(segTwo.getEnd()-segOne.getStart()) < 1000 and self.closeOnContig(segOne, segTwo)) or ((segTwo.getEnd() > segOne.getStart())and (segTwo.getStart() < segOne.getStart())):
            newSegment=copy.copy(segTwo)
            newSegment.relEnd=segOne.getRelEnd()
            newSegment.end=segOne.getEnd()
            return newSegment
        elif segOne.getStart() == segTwo.getStart() and segOne.getEnd() ==segTwo.getEnd():
            return segOne
        else:
            return False

    def closeOnContig(self, segOne, segTwo):
        if segOne.getStrand() == '+':
            if abs(segOne.getConEnd() - segTwo.getConStart()) <1000:
                return True
        elif segOne.getStrand() == '-':
            if abs(segTwo.getConEnd() - segOne.getConStart()) <1000:
                return True
        else:
            return False
    def getLengthInScaffold(self,scaffoldName):
        '''outputs an integer of how many bases of the contig map to the scaffold'''
        connectorNames=[s.getName() for s in self.getConnectors()]
        if scaffoldName not in connectorNames:
            return 0
        else:
            totalLength=0
            for seg in self.getGoodSegments():
                if seg.getOverlap()[0].getName() == scaffoldName:
                    totalLength+=seg.getLength()
            return totalLength
