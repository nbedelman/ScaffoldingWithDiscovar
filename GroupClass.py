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

class Group(object):
    def __init__(self, contigList, scaffoldList):
        self.contigList=contigList
        self.scaffoldList=scaffoldList
        self.superScaffolds=[]
        self.fullGroupJoin=[]
    def getContigList(self):
        return self.contigList
    def getScaffoldList(self):
        return self.scaffoldList
    def getSuperScaffolds(self):
        return self.superScaffolds
    def getFullGroupJoin(self):
        return self.fullGroupJoin
    def printFullGroupJoin(self):
        output=''
        for part in self.getFullGroupJoin().getPartsInOrder():
            output+=part.printPart()
        return output
    def printGroupEnvelopers(self):
        output=''
        envelopeDict= self.getFullGroupJoin().getEnveloperDict()
        for key in envelopeDict:
            ground=key
            allEnveloped=envelopeDict[key][0]
            enveloped=''
            for env in allEnveloped:
                enveloped+=str(env)
            pairs=''
            for scaffold in envelopeDict[key][1:]:
                pairs += "(%s, %s)" % (scaffold[0],scaffold[1])
            output+="%s : %s (%s)\n" % (ground, enveloped,pairs)
        return output
            
    def makeSuperScaffolds(self):
        superScaffolds=[]
        hiddenIgnores=[]
        for contig in self.getContigList():
            segments=[]
            smallIgnores=[]
            for seg in contig.getCombinedSegments():
                if seg.getLength() < 1000:
                    smallIgnores.append(copy.copy(seg))
                else:
                    segments.append(copy.copy(seg))
            orderedSegs=self.orderSegs(segments)
            palindromeChecked=self.checkPalindrome(orderedSegs)
            joinedSegments=self.joinEachSegment(palindromeChecked, [],smallIgnores)
            joinedSegs=joinedSegments[0]
            ignoredSegs=joinedSegments[1]
            if joinedSegs==[]:
                hiddenIgnores+=ignoredSegs
            fullContigJoin=self.joinSuperScaffolds(joinedSegs,[])
            if fullContigJoin != []:
                if fullContigJoin.getContigs() == []:
                    fullContigJoin.contigs.append(contig)
                self.superScaffolds.append(fullContigJoin)
                superScaffolds.append(fullContigJoin)
        toOrder=[(sup, sup.getTotalLength()) for sup in superScaffolds]
        ordered=sorted(toOrder, key=itemgetter(1), reverse=True)
        orderedSupers=[sup[0] for sup in ordered]
        fullGroupJoin=self.joinSuperScaffolds(orderedSupers,[])
        self.fullGroupJoin=fullGroupJoin
        if self.fullGroupJoin != []:
            usedScafs=[scaf.getName() for scaf in self.fullGroupJoin.getUsedScaffolds()]
            unusedScafs=[]
            for scaf in self.getScaffoldList():
                if scaf.getName() not in usedScafs:
                    unusedScafs.append(scaf)
            self.fullGroupJoin.alignWithBestScaf()
            self.fullGroupJoin.segIgnores+=hiddenIgnores
            self.fullGroupJoin.unusedScaffolds+=unusedScafs
        
        
    def joinSuperScaffolds(self, superScaffolds, joinedSupers):
        if joinedSupers == [] and len(superScaffolds) > 0:
            joinedSupers=superScaffolds[0]
            return self.joinSuperScaffolds(superScaffolds[1:],joinedSupers)
        elif len(superScaffolds) == 0:
            return joinedSupers
        else:
            notJoined=[]
            for index in range(len(superScaffolds)):
                toJoin=superScaffolds[index]
                if toJoin.isOverlapping(joinedSupers):
                    if toJoin.hasNewInfo(joinedSupers):
                        joinedSuperScafs=joinedSupers.getUsedScaffolds()
                        joinedSuperContigs=joinedSupers.getContigs()
                        
                        overlappingPart=toJoin.getFirstOverlap(joinedSupers)
                        toJoin.makePositive(overlappingPart)
                        joinedSupers.makePositive(overlappingPart)
                        supers=[joinedSupers,toJoin]
                        
                        joinedSupersStart=joinedSupers.lengthBefore(overlappingPart)
                        toJoinStart=toJoin.lengthBefore(overlappingPart)
                        distFromStart=[joinedSupersStart,toJoinStart]
                        
                        joinedSupersEnd=joinedSupers.lengthAfter(overlappingPart)
                        toJoinEnd=toJoin.lengthAfter(overlappingPart)                        
                        distFromEnd=[joinedSupersEnd,toJoinEnd]
                        
                        furthestFromStart=distFromStart.index(max(distFromStart))
                        startSuper=supers[furthestFromStart]
                        
                        furthestFromEnd=distFromEnd.index(max(distFromEnd))
                        endSuper=supers[furthestFromEnd]
                        
                        startOverlappingIndex=startSuper.getOverlappingIndices(overlappingPart)[0]
                        endOverlappingIndex=endSuper.getOverlappingIndices(overlappingPart)[0]
                        
                        newStart=[]
                        for part in range(len(startSuper.getPartsInOrder()[:startOverlappingIndex])):
                            newStart.append(startSuper.getPartsInOrder()[part].exportPart())
                            
                        newEnd=[]
                        try:
                            for part in range(endOverlappingIndex+1, len(endSuper.getPartsInOrder())):
                                newEnd.append(endSuper.getPartsInOrder()[part].exportPart())
                        except IndexError:
                            newEnd=[]
                        
                        
                        startOverlapPart=startSuper.getPartsInOrder()[startOverlappingIndex]
                        endOverlapPart=endSuper.getPartsInOrder()[endOverlappingIndex]
                        overlapBackbone=startOverlapPart.getBackbone()
                        
                        newOverlapPart=[(overlapBackbone,max(startOverlapPart.getStart(), endOverlapPart.getStart()), min(startOverlapPart.getEnd(), endOverlapPart.getEnd()),'+'),]
                        joinedSupers=newStart+newOverlapPart+newEnd
                        joinedSupers=self.checkNegatives(joinedSupers)
                        joinedSupers = SuperSegment(joinedSupers)
                        if toJoin.getContigs() != []:
                            joinedSupers.contigs.append(toJoin.contigs[0])
                        joinedSupers.contigs+=joinedSuperContigs
                        joinedSupers.usedScaffolds+=toJoin.getUsedScaffolds()
                        joinedSupers.usedScaffolds+=joinedSuperScafs
                        
                        try:
                            unjoined=notJoined+superScaffolds[index+1:]
                        except IndexError:
                            unjoined=notJoined
                        return self.joinSuperScaffolds(unjoined, joinedSupers)
                    else:
                        if toJoin.getContigs() != []:
                            joinedSupers.contigs.append(toJoin.contigs[0])
                        joinedSupers.usedScaffolds+=toJoin.getUsedScaffolds()
                        try:
                            unjoined=notJoined+superScaffolds[index+1:]
                        except IndexError:
                            unjoined=notJoined
                        return self.joinSuperScaffolds(unjoined, joinedSupers) 
                else:
                    notJoined.append(superScaffolds[index]) 
        #If there are superScaffolds that could not be joined, it's because we have something like
        #scaffolds 1 and 2 both were placed inside 3, and the superScaffold that links 1 and 2 together
        #therefore doesn't overlap with the big combined one
        for supScaf in notJoined:
            if supScaf.getContigs() != []:
                joinedSupers.contigs.append(supScaf.getContigs()[0])
            joinIgnores=[seg.getOverlap()[0].getName() for seg in joinedSupers.getTrueIgnores()]
            for used in supScaf.getUsedScaffolds():
                if used.getName() not in joinIgnores:
                    joinedSupers.usedScaffolds.append(used)
        return joinedSupers   
            
        
    def checkPalindrome(self,segmentList):
        tempInfo=[]
        dePalindromed=[]
        for segment in segmentList:
            tempInfo.append((segment,segment.getScafStart(), segment.getScafEnd(), segment.getStrand()))
        for info in range(len(tempInfo)-1):
            if tempInfo[info][1:2]==tempInfo[info+1][1:2]:
                if tempInfo[info][3]=='+':
                    dePalindromed.append(tempInfo[info][0])
                else:
                    dePalindromed.append(tempInfo[info+1][0])
        if len(dePalindromed)!=0:
            return dePalindromed
        else:
            return segmentList
            
    def orderSegs(self,segmentList):
        unordered=[]
        ordered=[]
        for i in segmentList:
            unordered.append((i,i.getConStart()))
        ordered=sorted(unordered, key=itemgetter(1))
        onlySegs=[segment[0] for segment in ordered]
        return onlySegs
        
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
            result=self.joinSegments(joinThisTime)
            if result:
                updatedJoinedSegments.append(result)
        except NestedSegsError:
            ignoredSegs+=[segmentList[joinIndexTwo]]
            result=self.nestedSegBranch(segmentList,joinIndexTwo + 1, ignoredSegs)
            if result[0]:
                updatedJoinedSegments.append(result[0])
            ignoredSegs+=result[2]
            joinIndexTwo=result[1]   
        return self.joinEachSegment(segmentList[joinIndexTwo:], updatedJoinedSegments, ignoredSegs)
            
    def nestedSegBranch(self,segmentList, indTwo, ignoredSegs):
        #print "nestedSegBranch"
        try:
            joinThisTime=[segmentList[0],segmentList[indTwo]]
            output=self.joinSegments(joinThisTime)
            return (output, indTwo, ignoredSegs)
        except NestedSegsError:
            ignoredSegs.append(segmentList[indTwo])
            return self.nestedSegBranch(segmentList, indTwo+1, ignoredSegs)
        except IndexError:
            return (None, indTwo, ignoredSegs)
    
    
    def joinSegments(self,segmentPair):
        seg1=copy.copy(segmentPair[0])
        seg2=copy.copy(segmentPair[1])
        if seg1.getStrand()==seg2.getStrand():
            return self.sameDirectionBranch(seg1,seg2)
        else:
            return self.diffDirectionBranch(seg1,seg2)
    
    def sameDirectionBranch(self, seg1,seg2):
        #print "sameDirectionBranch"
        if seg1.getOverlap()[0].getName() == seg2.getOverlap()[0].getName():
            output=[(seg1.getOverlap()[0], 1, seg1.getOverlap()[0].getLength(), seg1.getOverlap()[0].getStrand()),]
            output=self.checkNegatives(output)
            output=SuperSegment(output)
            output.contigs.append(seg1.getContig())
            output.usedScaffolds+=[seg1.getOverlap()[0],]
            return output
                       
        else:
            return self.differentScaffoldBranch(seg1,seg2)
    
    def diffDirectionBranch(self,seg1,seg2):
        #print "diffDirectionBranch"
        if seg1.getOverlap()[0].getScore() > seg2.getOverlap()[0].getScore():
            newSegTwo = copy.copy(seg2)
            newSegTwo.flipSegment()
            return self.sameDirectionBranch(seg1, newSegTwo)
        elif (seg1.getOverlap()[0].getScore() == seg2.getOverlap()[0].getScore()) and (seg1.getOverlap()[0].getLength() > seg2.getOverlap()[0].getLength()):
            newSegTwo = copy.copy(seg2)
            newSegTwo.flipSegment()
            return self.sameDirectionBranch(seg1, newSegTwo)
        else:
            newSegOne = copy.copy(seg1)
            newSegOne.flipSegment()
            return self.sameDirectionBranch(newSegOne, seg2)
        
        
        
    def differentScaffoldBranch(self,seg1,seg2):
        #print "differentScaffoldBranch"
        if seg1.getConEnd() < seg2.getConEnd():
            return self.diffScaffCombine(seg1,seg2) #was notNestedSegsBranch(seg1,seg2)
        else:
            #print "nestedSegs"
            raise NestedSegsError('''Contig %s has nested segments %s and %s''' % (seg1.getContig().getName(), seg1.getName(), seg2.getName()))

    def diffScaffCombine(self,seg1,seg2):
        distFromScafStart=[]
        distFromScafEnd=[]
        scaffolds=[copy.copy(seg1.getOverlap()[0]),copy.copy(seg2.getOverlap()[0])]
        contig=copy.copy(seg1.getContig())
        output=''
        if seg1.getStrand() == '+' and self.distanceBetweenSegsChecks(seg1,seg2):
            distFromScafStart.append(seg1.getDistanceFromScafStart())
            distFromScafStart.append(seg2.getDistanceFromScafStart() - (seg2.getConStart()-seg1.getConStart()))
            distFromScafEnd.append(seg1.getDistanceFromScafEnd() - seg2.getConEnd()-seg1.getConEnd())
            distFromScafEnd.append(seg2.getDistanceFromScafEnd())
            furthestFromStart=distFromScafStart.index(max(distFromScafStart))
            startScaf=scaffolds[furthestFromStart]
            furthestFromEnd=distFromScafEnd.index(max(distFromScafEnd))
            endScaf=scaffolds[furthestFromEnd]
            
        elif seg1.getStrand() == '-' and self.distanceBetweenSegsChecks(seg1,seg2):
            distFromScafStart.append(seg1.getDistanceFromScafStart()- (seg2.getConEnd()-seg1.getConEnd()))
            distFromScafStart.append(seg2.getDistanceFromScafStart())
            distFromScafEnd.append(seg1.getDistanceFromScafEnd())
            distFromScafEnd.append(seg2.getDistanceFromScafEnd()- seg2.getConStart()-seg1.getConStart())
            furthestFromStart=distFromScafStart.index(max(distFromScafStart))
            startScaf=scaffolds[furthestFromStart]
            furthestFromEnd=distFromScafEnd.index(max(distFromScafEnd))
            endScaf=scaffolds[furthestFromEnd]

        else:
            return  
        
        if startScaf.getName() == endScaf.getName():
                output=[(startScaf, 1, startScaf.getLength(), startScaf.getStrand()),]
                output=self.checkNegatives(output)
                output=SuperSegment(output)
                output.contigs.append(contig)
                output.usedScaffolds+=[seg1.getOverlap()[0],seg2.getOverlap()[0]]
                
        else:
            endOfStartScafOverlap = min(startScaf.getLength()-distFromScafEnd[furthestFromStart], startScaf.getLength())
            firstScafSeg=(startScaf,1,endOfStartScafOverlap, startScaf.getStrand())
            
            contigDiff=seg2.getConStart()-seg1.getConEnd()
            if contigDiff <0:
                startOfEndScafOverlap=max(distFromScafStart[furthestFromEnd] + contigDiff, 1)
                endScafSeg=(endScaf,startOfEndScafOverlap, endScaf.getLength(), endScaf.getStrand())
                output=[firstScafSeg,endScafSeg]
            elif contigDiff >=0:
                startOfEndScafOverlap=max(distFromScafStart[furthestFromEnd], 1)
                contigSeg=(contig,min(seg1.getConEndPos(), seg2.getConEndPos()), max(seg1.getConStartPos(), seg2.getConStartPos()),seg1.getStrand())
                endScafSeg=(endScaf, startOfEndScafOverlap, endScaf.getLength(), endScaf.getStrand())
                output=[firstScafSeg,contigSeg,endScafSeg]                              
            output=self.checkNegatives(output)
            output=SuperSegment(output)
            output.usedScaffolds+=[seg1.getOverlap()[0],seg2.getOverlap()[0]]
        return output
                    
    
    def distanceBetweenSegsChecks(self,seg1, seg2):
        #print "distanceBetweenSegsChecks"
        return seg1.getConEnd()-seg2.getConEnd() < 1000
                    

    
    def checkNegatives(self, joinedSeg):
        newSeg=copy.copy(joinedSeg)
        for entry in range(len(newSeg)):
            newSeg[entry]=list(newSeg[entry])
            for number in range(len(newSeg[entry])):
                if newSeg[entry][number] < 1:
                    newSeg[entry][number] = 1
            newSeg[entry]=tuple(newSeg[entry])
        return newSeg
                    

            
           