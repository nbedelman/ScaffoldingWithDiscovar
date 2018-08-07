#InversionCandidateClass
from Bio.Seq import Seq
from Bio import SeqIO
import os
import itertools
import string
import copy
from operator import itemgetter
from ScaffoldClass import *
from ContigClass import *
from SegmentClass import *
from SuperSegmentClass import *
from PartClass import *
from ErrorClasses import *

class InversionCandidate(Contig):
    def __init__(self, contig):
        self.contig=contig
        self.isInversion=False
        self.positiveSegments=[]
        self.negativeSegments=[]
        self.mainStrand=self.contig.getStrand()
        self.segments=self.contig.getCombinedSegments()

    def getContig(self):
        return self.contig

    def getSegments(self):
        return self.segments
    def getMainStrand(self):
        return self.mainStrand

    def findInversion(self, pctThreshold=.1, minLength=1000):
        segments=[]
        smallIgnores=[]
        for seg in self.getSegments():
            if seg.getLength() < minLength:
                smallIgnores.append(copy.copy(seg))
            else:
                segments.append(copy.copy(seg))
        orderedSegs=self.orderSegs(segments)
        for seg in orderedSegs:
            if seg.getStrand()!=self.getMainStrand():
                sizeAndContiguityCheck=self.verifySizeAndContiguity(orderedSegs, pctThreshold)
                overlapCheck=self.verifyNoOverlapsAndSameScaf(orderedSegs)
                break
        try:
            self.isInversion = sizeAndContiguityCheck and overlapCheck
        except UnboundLocalError:
            self.isInversion = False

    def verifyNoOverlapsAndSameScaf(self,orderedSegs):
        for i in itertools.combinations(orderedSegs,2):
            if self.areOverlapping(i[0],i[1]):
                return False
            try:
                if i[0].getOverlap()[0] != i[1].getOverlap()[0]:
                    return False
            except IndexError:
                pass
        return True

    def verifySizeAndContiguity(self, segList, pctThreshold):
        ''' if we have an indication of inversion, check to make sure at least 10% of the contig is aligning in each direction
        In addition, there should be at most one 'change of direction', which would be if this was in internal inversion.
        Make sure that is the case.'''
        numChanges=-1
        currentDirection=''
        positiveLength=0
        negativeLength=0
        internalOK=True
        for seg in segList:
            if seg.getStrand()=='+':
                positiveLength+=seg.getLength()
            else:
                negativeLength+=seg.getLength()
            if seg.getStrand()!=currentDirection:
                currentDirection=seg.getStrand()
                numChanges+=1
        posOK=float(positiveLength)/self.getContig().getLength() > pctThreshold
        negOK=float(negativeLength)/self.getContig().getLength() > pctThreshold
        fullOK=(float(positiveLength)/self.getContig().getLength() + float(negativeLength)/self.getContig().getLength()) > .5
        if numChanges==2:
            internalOK=self.isTrueInternal(segList)
        return posOK and negOK and internalOK and fullOK and (numChanges < 3)

    def isTrueInternal(self,segList):
        goodInternal=True
        mainDirection=segList[0].getStrand()
        scafStart=segList[0].getRelStart()
        scafEnd=segList[0].getRelEnd()
        for seg in segList:
            if (seg.getStart()<scafStart) or (seg.getEnd()>scafEnd):
                goodInternal=False
        return goodInternal

    def areOverlapping(self,seg1,seg2):
        if (seg1.getConStart() >= seg2.getConStart()) and (seg1.getConStart() <= seg2.getConEnd()):
            return True
        elif (seg2.getConStart() >= seg1.getConStart()) and (seg2.getConStart() <= seg1.getConEnd()):
            return True
        else:
            return False

    def getBreakPoints(self):
        segList=[]
        segments=self.getSegments()
        for seg in segments:
            segList.append([seg.getConStart(), seg.getConEnd(), seg.getOriginalStrand(),seg.getStart(), seg.getEnd()])
        segList=sorted(segList, key=itemgetter(0))
        startDir=segList[0][2]
        for entry,seg in enumerate(segList):
            if seg[2] != startDir:
                firstSeg=segList[entry-1]
                secondSeg=segList[entry]
                break
        marginOfError=secondSeg[0]-firstSeg[1]
        if startDir=='+':
            bp1Start=firstSeg[4]
            bp1End=bp1Start+marginOfError
            bp2Start=secondSeg[4]
            bp2End=bp2Start+marginOfError
        else:
            bp1End=firstSeg[3]
            bp1Start=bp1End-marginOfError
            bp2End=secondSeg[3]
            bp2Start=bp2End-marginOfError
        bp1='''%s\t%i\t%i\t%s\n''' % (self.getContig().getChrom(),bp1Start,bp1End,self.getContig().getName())
        bp2='''%s\t%i\t%i\t%s\n''' % (self.getContig().getChrom(),bp2Start,bp2End,self.getContig().getName())
        return (bp1+bp2)
