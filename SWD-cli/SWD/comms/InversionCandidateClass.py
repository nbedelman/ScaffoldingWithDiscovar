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
    
    def findInversion(self):
        segments=[]
        smallIgnores=[]
        for seg in self.getSegments():
            if seg.getLength() < 1000:
                smallIgnores.append(copy.copy(seg))
            else:
                segments.append(copy.copy(seg))
        orderedSegs=self.orderSegs(segments)
        for seg in orderedSegs:
            if seg.getStrand()!=self.getMainStrand():
                sizeAndContiguityCheck=self.verifySizeAndContiguity(orderedSegs)
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

    def verifySizeAndContiguity(self, segList):
        ''' if we have an indication of inversion, check to make sure at least 10% of the contig is aligning in each direction
        In addition, there should be at most one 'change of direction', which would be if this was in internal inversion. 
        Make sure that is the case.'''
        numChanges=-1
        currentDirection=''
        positiveLength=0
        negativeLength=0
        for seg in segList:
            if seg.getStrand()=='+':
                positiveLength+=seg.getLength()
            else:
                negativeLength+=seg.getLength()
            if seg.getStrand()!=currentDirection:
                currentDirection=seg.getStrand()
                numChanges+=1
        posOK=float(positiveLength)/self.getContig().getLength() > .1
        negOK=float(negativeLength)/self.getContig().getLength() > .1
        return posOK and negOK and numChanges<3
    
    def areOverlapping(self,seg1,seg2):
        if (seg1.getConStart() > seg2.getConStart()) and (seg1.getConStart() < seg2.getConEnd()):
            return True
        elif (seg2.getConStart() > seg1.getConStart()) and (seg2.getConStart() < seg1.getConEnd()):
            return True
        else:
            return False
        
    
    def joinScaffolds(self):
        segments=[]
        smallIgnores=[]
        hiddenIgnores=[]
        for seg in self.getContig().getCombinedSegments():
            if seg.getLength() < 1000:
                smallIgnores.append(copy.copy(seg))
            else:
                segments.append(copy.copy(seg))
        orderedSegs=self.orderSegs(segments)
        try:
            joinedSegments=self.joinEachSegment(orderedSegs, [],smallIgnores)
            joinedSegs=joinedSegments[0]
            ignoredSegs=joinedSegments[1]
            if joinedSegs==[]:
                hiddenIgnores+=ignoredSegs
            fullJoin=self.joinSuperScaffolds(joinedSegs,[], intraContig=True)
            return fullJoin, hiddenIgnores
        except InversionError:
            raise
            
    def joinEachSegment(self, segmentList, joinedSegments, ignoredSegs):
        if len(segmentList)<=1:
            if  joinedSegments != []:
                for joined in joinedSegments:
                    joined.segIgnores+=ignoredSegs
            return (joinedSegments,ignoredSegs)
        else:
            updatedJoinedSegments=joinedSegments
            joinIndexOne=0
            joinIndexTwo=1

        try:
            joinThisTime=[segmentList[joinIndexOne],segmentList[joinIndexTwo]]
            result=self.joinSegments(joinThisTime, "exclusive")
            if result:
                result.inclusiveJoin=self.joinSegments(joinThisTime, "inclusive")
                updatedJoinedSegments.append(result)
        except NestedSegsError:
            ignoredSegs+=[segmentList[joinIndexTwo]]
            result=self.nestedSegBranch(segmentList,joinIndexTwo + 1, ignoredSegs, "exclusive")
            if result[0]:
                result[0].inclusiveJoin=self.nestedSegBranch(segmentList,joinIndexTwo + 1, ignoredSegs, "inclusive")[0]
                updatedJoinedSegments.append(result[0])
            ignoredSegs+=result[2]
            joinIndexTwo=result[1]
        except InversionError:
            raise
        return self.joinEachSegment(segmentList[joinIndexTwo:], updatedJoinedSegments, ignoredSegs)

    def orderSegs(self,segmentList):
        unordered=[]
        ordered=[]
        for i in segmentList:
            unordered.append((i,i.getConStart()))
        ordered=sorted(unordered, key=itemgetter(1))
        onlySegs=[segment[0] for segment in ordered]
        return onlySegs
    
    def joinSuperScaffolds(self, superScaffolds, joinedSupers, intraContig):
        if joinedSupers == [] and len(superScaffolds) > 0:
            joinedSupers=superScaffolds[0]
            return self.joinSuperScaffolds(superScaffolds[1:],joinedSupers,intraContig)
        elif len(superScaffolds) == 0:
            return joinedSupers
        else:
            notJoined=[]
            for index in range(len(superScaffolds)):
                toJoin=superScaffolds[index]
                if toJoin.isOverlapping(joinedSupers):
                    try:
                        together=self.overlapJoin(toJoin, joinedSupers, intraContig)
                        try:
                            unjoined=notJoined+superScaffolds[index+1:]
                        except IndexError:
                            unjoined=notJoined
                        return self.joinSuperScaffolds(unjoined, together, intraContig)
                    except UnboundLocalError:
                        #This happens in some complex enveloper cases. There are no shared scaffolds between the superScaffolds, but the contig maps to both of them. Here, I essentially save the join for the next iteration.
                        if joinedSupers.getFirstOverlap(toJoin, partType=Contig):
                            notJoined.append(superScaffolds[index])
                    except AmbiguousJoinError:
                        #Another complex enveloper case; I don't want to include the second superScaffold(here superScaffolds[index]), so I can't put it in notJoined like before. Just skip it.
                        pass
                else:
                    notJoined.append(superScaffolds[index])

            if notJoined != []:
                split=self.insideOrUnjoinable(joinedSupers, notJoined)
                joinedSupers=split[0]
                unJoinable=split[1]
                if (unJoinable != []) and (not intraContig):
                    #only want to make new groups if we're combining contigs. Each contig is only allowed to to one thing per round!
                    allContigs=[]
                    allScaffolds=[]
                    for unJoined in unJoinable:
                        allContigs+=unJoined.getContigs()
                        allScaffolds+=unJoined.getUsedScaffolds()
                    newGroup = Group(allContigs, allScaffolds, self.chromosome)
                    newGroup.brandNew=False
                    newGroup.transferredSuperScaffolds=unJoinable
                    self.chromosome.groups.append(newGroup)
            return joinedSupers

    def nestedSegBranch(self,segmentList, indTwo, ignoredSegs, joinType):
        #print "nestedSegBranch"
        try:
            joinThisTime=[segmentList[0],segmentList[indTwo]]
            output=self.joinSegments(joinThisTime, joinType)
            return (output, indTwo, ignoredSegs)
        except NestedSegsError:
            ignoredSegs.append(segmentList[indTwo])
            return self.nestedSegBranch(segmentList, indTwo+1, ignoredSegs, joinType)
        except IndexError:
            return (None, indTwo, ignoredSegs)


    def joinSegments(self,segmentPair, joinType):
        seg1=copy.copy(segmentPair[0])
        seg2=copy.copy(segmentPair[1])
        if seg1.getOriginalStrand()==seg2.getOriginalStrand():
            return self.sameDirectionBranch(seg1,seg2, joinType)
        else:
            return self.diffDirectionBranch(seg1,seg2, joinType)


    def sameDirectionBranch(self, seg1,seg2, joinType):
        #print "sameDirectionBranch"
        if seg1.getOverlap()== seg2.getOverlap()[0]:
            output=[(seg1.getOverlap()[0], 0, seg1.getOverlap()[0].getLength(), seg1.getOverlap()[0].getStrand()),]
            output=self.checkNegatives(output)
            output=SuperSegment(output)
            output.contigs.append(seg1.getContig())
            output.usedScaffolds+=[seg1.getOverlap()[0],]
            return output
        else:
            return False

    def diffDirectionBranch(self,seg1,seg2, joinType):
        #print "diffDirectionBranch"
        if seg1.getOverlap()[0] == seg2.getOverlap()[0]:
            self.isInversion=True
            return False
            #raise InversionError('''contig %s has a small inversion relative to the reference. please adjust manually.''' % (seg1.getContig().getName()))
        

    def distanceBetweenSegsChecks(self,seg1, seg2):
        #print "distanceBetweenSegsChecks"
        #See if segments are next to each other on the contig that they come from.
        #updated June 26, 2017. Not really sure what a reasonable number is for this.
        return abs(seg2.getConStart()-seg1.getConEnd()) < 10000



    def checkNegatives(self, joinedSeg):
        newSeg=copy.copy(joinedSeg)
        for entry in range(len(newSeg)):
            newSeg[entry]=list(newSeg[entry])
            for number in range(len(newSeg[entry])):
                if newSeg[entry][number] < 0:
                    newSeg[entry][number] = 0
            newSeg[entry]=tuple(newSeg[entry])
        return newSeg
