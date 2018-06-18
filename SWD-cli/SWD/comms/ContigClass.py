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

    def combineSegments(self,multiScafs=True):
        '''here, combine segments that are very close to each other and really represent a single long mapping.
        Also get rid of very small (and therefore possibly spurious) alignments. Here, "small" is <500 bp AND < 75% of the scaffold to which it maps'''
        if self.getGoodSegments() == []:
            raise NoSegmentError("Remember to cull segments before combining!")
        #elif len(self.getGoodSegments()) ==1:
        #    raise NotInformativeError("After culling, this contig no longer maps to two scaffolds")
        else:
            orderedSegs=self.orderSegs(self.getGoodSegments())
            combinedSegs=self.combineLists(orderedSegs, multiScafs)
            #if len(combinedSegs) == 1:
            #    raise NotInformativeError("After combining, this contig no longer maps to two scaffolds")
            #else:
            longEnough=[]
            for seg in combinedSegs:
                #Added the multiScafs flag when creating the findInversions functionality. 
                #In findInversions, we're using sincle-contig data only.
                if multiScafs:
                    overlapLength=seg.getOverlap()[0].getLength()
                    mappingPercent=float(seg.getConEnd()-seg.getConStart())/overlapLength
                else:
                    mappingPercent=0
                if (seg.getConEnd()-seg.getConStart() > 500) or (mappingPercent > .75) :
                    longEnough.append(seg)
            self.combinedSegments=longEnough
    def outputBed(self, segList, fileName):
        with open(fileName, "a+") as o:
            for segment in segList:
                outList=[]
                outString=''
                for i in [segment.getChrom(),segment.getRelStart(),segment.getRelEnd(), segment.getName(), segment.getScore()\
                ,segment.getStrand(), segment.getStart(), segment.getEnd()]:
                    outString+=str(i)+"\t"
                outString+=segment.getColor()+"\n"
                o.write(outString)
        o.close()
    def findConnectors(self, scaffoldDict, segType):
        connects=[]
        removeIndices=[]
        if segType=='good':
            segList=self.getGoodSegments()
            markScaffolds=False
        elif segType=='combined':
            segList=self.getCombinedSegments()
            markScaffolds=True
        for idx, segment in enumerate(segList):
            #If the segment maps uniquely, this will return the connecting scaffold.
            #If it does not, this will return the segment split into two at the scaffold break
            #These have to be run through the process again.
            overlapper=segment.findOverlaps(scaffoldDict, markScaffolds)
            if len(overlapper)==1:
                scaffold=overlapper[0]
                if (not scaffold in connects) and (not scaffold.getName() == "Ns"):
                    connects.append(scaffold)
            else:
                #This can only happen the first time, which means the segList will be goodSegments.
                segList+=overlapper
                #self.goodSegments+=overlapper
                removeIndices.append(idx)
        
        #If there are multiple to remove, must remove the largest first so it doesn't screw up the order.
        for index in sorted(removeIndices, reverse=True):
            del self.goodSegments[index]
        
        self.connectors=connects

    def orderSegs(self,segmentList):
        unordered=[]
        ordered=[]
        for i in segmentList:
            unordered.append((i,i.getStart()))
        ordered=sorted(unordered, key=itemgetter(1))
        segsOnly=[entry[0] for entry in ordered]
        return segsOnly



    def combineLists(self, goodSegs, multiScafs):
        combined=[]
        goodIndices=[]
        for i in (range(len(goodSegs)-1)):
            j=goodSegs[i]
            k=goodSegs[i+1]
            out = self.tryCombining(j,k, multiScafs)
            if out:
                goodIndices+=i,i+1
                combined.append(out)
        for l in range(len(goodSegs)):
            if not l in goodIndices:
                combined.append(goodSegs[l])
        if len(goodSegs) == len(combined):
            return goodSegs
        else:
            return self.combineLists(combined,multiScafs)


    def tryCombining(self, segOne, segTwo, multiScafs):
        if segOne.getStrand()!=segTwo.getStrand():
            return False
        elif multiScafs and (segOne.getOverlap()[0].getName() != segTwo.getOverlap()[0].getName()):
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
