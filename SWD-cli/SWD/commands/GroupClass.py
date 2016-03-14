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
    def __init__(self, contigList, scaffoldList, chromosome):
        self.contigList=contigList
        self.scaffoldList=scaffoldList
        self.superScaffolds=[]
        self.fullGroupJoin=[]
        self.chromosome=chromosome
        self.inversions=[]
        self.contigReports=[]
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
                enveloped+= '''%s, ''' % (str(env))
            pairs=''
            for scaffold in envelopeDict[key][1:]:
                pairs += "(%s, %s)" % (scaffold[0],scaffold[1])
            output+="%s : %s (%s)\n" % (ground, enveloped,pairs)
        return output
    def getContigReports(self):
        return self.contigReports        
    def makeSuperScaffolds(self):
        superScaffolds=[]
        hiddenIgnores=[]
        for contig in self.getContigList():
            try:
                fullContigJoin = self.joinFullContig(contig)
                contigJoin = fullContigJoin[0]
                contigIgnores = fullContigJoin [1]
                hiddenIgnores+=contigIgnores
                if contigJoin != []:
                    if contigJoin.getContigs() == []:
                        contigJoin.contigs.append(contig)
                    self.superScaffolds.append(contigJoin)
                    superScaffolds.append(contigJoin)
                    contigReport=self.reportContig(contig,contigJoin)
                    self.contigReports.append(contigReport)
            except InversionError:
                self.inversions.append(contig)
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
    
    def reportContig(self,contig,contigJoin):
        #make a report of what happened with each contig. List the contig, the scaffolds involved, and the final result. 
        #Indicate whether there was a join, an enveloping, an inversion, or nothing
        name=contig.getName()
        connectingScaffolds=[con.getName() for con in contig.getConnectors()]
        numConnects=len(connectingScaffolds)
        usedScaffolds=contigJoin.getUsedScaffolds()
        numUsed=len(usedScaffolds)
        finalSup=contigJoin.printSuperSeg()
        scafsInFinal=[]
        for part in contigJoin.getPartsInOrder():
            if part.getType()==Scaffold:
                scafsInFinal.append(part.getBackbone().getName())
        numScafsInFinal=len(scafsInFinal)
        if numUsed==numScafsInFinal:
            result="j"
        elif numScafsInFinal==1:
            result="e"
        else:
            result="ej"
        enveloped=""
        envelopeLength=""
        if "e" in result:
            for i in usedScaffolds:
                if i.getName() not in scafsInFinal:
                    enveloped+=(str(i.getName())+",")
                    envelopeLength+=(str(i.getLength())+",")
        lengths=[scaf.getLength() for scaf in contigJoin.getUsedScaffolds()]
        return [name,[scaf.getName() for scaf in usedScaffolds],lengths,enveloped,envelopeLength,result,"use"]
                
                
        
        
    def joinFullContig(self,contig):
        segments=[]
        smallIgnores=[]
        hiddenIgnores=[]
        for seg in contig.getCombinedSegments():
            if seg.getLength() < 1000:
                smallIgnores.append(copy.copy(seg))
            else:
                segments.append(copy.copy(seg))
        orderedSegs=self.orderSegs(segments)
        palindromeChecked=self.checkPalindrome(orderedSegs)
        if palindromeChecked != orderedSegs:
            hiddenIgnores+=palindromeChecked
            return [], hiddenIgnores
        else:
            try:
                joinedSegments=self.joinEachSegment(palindromeChecked, [],smallIgnores)
                joinedSegs=joinedSegments[0]
                ignoredSegs=joinedSegments[1]
                if joinedSegs==[]:
                    hiddenIgnores+=ignoredSegs
                if contig.getInOrEx() == "e":
                    fullJoin=self.joinSuperScaffolds(joinedSegs,[])
                elif contig.getInOrEx() == "i":
                    incJoins=[seg.getInclusiveJoin() for seg in joinedSegs]
                    fullJoin=self.joinSuperScaffolds(incJoins,[])
                return fullJoin, hiddenIgnores
            except InversionError:
                raise
        
        
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
                    together=self.overlapJoin(toJoin, joinedSupers)
                    try:
                        unjoined=notJoined+superScaffolds[index+1:]
                    except IndexError:
                        unjoined=notJoined
                    return self.joinSuperScaffolds(unjoined, together)   
                else:
                    notJoined.append(superScaffolds[index]) 
                
            if notJoined != [] :
                split=self.insideOrUnjoinable(joinedSupers, notJoined)
                joinedSupers=split[0]
                unJoinable=split[1]
                if unJoinable != []:
                    allContigs=[]
                    allScaffolds=[]
                    for unJoined in unJoinable:
                        allContigs+=unJoined.getContigs()
                        allScaffolds+=unJoined.getUsedScaffolds()
                    newGroup = Group(allContigs, allScaffolds, self.chromosome)
                    self.chromosome.groups.append(newGroup)
                    for contig in allContigs:
                        for groupCon in self.getContigList():
                            if contig.getName() == groupCon.getName():
                                self.contigList.remove(groupCon)
                    for scaffold in allScaffolds:
                        for groupScaf in self.getScaffoldList():
                            if scaffold.getName() == groupScaf.getName():
                                self.scaffoldList.remove(groupScaf)
            return joinedSupers  
                 
    
    def overlapJoin(self, joinedSupers,toJoin):
        '''assuming joinedSupers and toJoin overlap each other'''                  
        if toJoin.hasNewInfo(joinedSupers):
            newInf= self.newInfoJoin(joinedSupers,toJoin)
        else:
            if toJoin.getContigs() != []:
                joinedSupers.contigs.append(toJoin.contigs[0])
            joinedSupers.usedScaffolds+=toJoin.getUsedScaffolds()
            return joinedSupers 
            
    def newInfoJoin(joinedSupers, toJoin):  
        '''assuming toJoin has some new information to add to joinedSupers''' 
        envelopedScafFlag=False
        allJoined=[]
        allToJoin=[]
        doubleJoined=[]
        doubleToJoin=[]
        for part in joinedSupers.getPartsInOrder():
            if part.getName() in allJoined:
                doubleJoined.append(part.getName())
            allJoined.append(part.getName())
        for part in toJoin.getPartsInOrder(): 
            if part.getName() in allToJoin:
                doubletoJoin.append(part.getName())
            allToJoin.append(part.getName())
        enveloper=''
        for double in doubleJoined:
            if double in allToJoin:
                envelopedScafFlag=True
                enveloper=double
        if not envelopedScafFlag:
            for double in doubleToJoin:
                if double in allJoined:
                    envelopedScafFlag=True
                    enveloper=double                     
        if envelopedScafFlag:
            return self.envelopedScafJoin(joinedSupers, toJoin)
        else:
            return self.notEnvelopedScafJoin(joinedSupers, toJoin, enveloper)
                    
    def notEnvelopedScafJoin(self, joinedSupers, toJoin):
        joinedSuperScafs=joinedSupers.getUsedScaffolds()
        joinedSuperContigs=joinedSupers.getContigs()
        
        supers=self.alignDirections(joinedSupers, toJoin)
        
        overlappingPart=joinedSupers.getFirstOverlap(toJoin, partType=Scaffold)
        distFromStart=[joinedSupers.lengthBefore(overlappingPart),toJoin.lengthBefore(overlappingPart)]
        
    
        distFromEnd=[joinedSupers.lengthAfter(overlappingPart),toJoin.lengthAfter(overlappingPart)]
        
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
        return joinedSupers
    
    def envelopedScafJoin(self, joinedSupers, toJoin, enveloper):
        #In the simplest case, one of the scaffolds is simply the entire scaffold (S1), and the other is S1, S2, S1
        if len(joinedSupers.getPartsInOrder()) ==1:
            if joinedSupers.getContigs() != []:
                toJoin.contigs.append(joinedSupers.contigs[0])
            toJoin.usedScaffolds+=joinedSupers.getUsedScaffolds()
            joinedSupers=toJoin
            return joinedSupers
        elif len(toJoin.getPartsInOrder())==1:
            if toJoin.getContigs() != []:
                joinedSupers.contigs.append(toJoin.contigs[0])
            joinedSupers.usedScaffolds+=toJoin.getUsedScaffolds()
            return joinedSupers
        else:
            return self.complexEnvScafJoin(joinedSupers, toJoin, enveloper)
    
    def complexEnvScafJoin(self, joinedSupers, toJoin, enveloper):
        #In the complex case, one or both of the scaffolds uses S1 as an enveloper, and both have multiple scaffolds
        #first, figure out where in each superScaffold the enveloper is
        supers=self.alignDirections(joinedSupers, toJoin)
        joinedSupersEnvIndices=[]
        for i in range(len(joinedSupers.getPartsInOrder())):
            if joinedSupers.getPartsInOrder()[i].getName() == enveloper:
                joinedSupersEnvIndices.append(i)
        
        toJoinEnvIndices=[]
        for j in range(len(toJoin.getPartsInOrder())):
            if toJoin.getPartsInOrder()[j].getName() == enveloper:
                toJoinEnvIndices.append(j)
        #next, as for the non-enveloping case, figure out which superScaffold goes first and which goes second (could be the same one)
        distFromStart=[joinedSupers.lengthBefore(min(joinedSupersEnvIndices)),toJoin.lengthBefore(min(toJoinEnvIndices))]  
        distFromEnd=[joinedSupers.lengthAfter(max(joinedSupersEnvIndices)),toJoin.lengthAfter(max(toJoinEnvIndices))]   
        furthestFromStart=distFromStart.index(max(distFromStart))
        startSuper=supers[furthestFromStart]
        furthestFromEnd=distFromEnd.index(max(distFromEnd))
        endSuper=supers[furthestFromEnd]
        #Finally, figure out                 
                                                            
    def insideOrUnjoinable(self,joinedSupers, notJoinedList):
        alreadyJoined=copy.copy(joinedSupers)
        candidateList=copy.copy(notJoinedList)
        startNum=len(candidateList)
        stillUnjoined=[]
        
        alreadyJoined.findEnvelopers()
        allInvolved=alreadyJoined.getEnveloped()+[used.getName() for used in alreadyJoined.getUsedScaffolds()]
        
        for supScaf in candidateList:
            reallyAbsent=True
            supScaf.findEnvelopers()
            allInUnjoined=supScaf.getEnveloped()+[used.getName() for used in supScaf.getUsedScaffolds()]
            for outScaf in allInUnjoined:
                for inScaf in allInvolved:
                    if outScaf == inScaf:
                        reallyAbsent=False
            
            if not reallyAbsent:
                superScafs=[alreadyJoined,supScaf]
                sizes=[alreadyJoined.getTotalLength(), supScaf.getTotalLength()]
                bigger=superScafs[sizes.index(max(sizes))]
                smaller=superScafs[sizes.index(min(sizes))]
                if smaller.getContigs() != []:
                    bigger.contigs+=supScaf.getContigs()
                joinIgnores=[seg.getOverlap()[0].getName() for seg in bigger.getTrueIgnores()]
                for used in smaller.getUsedScaffolds():
                    if used.getName() not in joinIgnores:
                        bigger.usedScaffolds.append(used)
                        allInvolved=bigger.getEnveloped()+[used.getName() for used in bigger.getUsedScaffolds()]
                alreadyJoined=bigger
            else:
                stillUnjoined.append(supScaf)
        if startNum == len(stillUnjoined):
            return alreadyJoined, stillUnjoined
        else:
            return self.insideOrUnjoinable(alreadyJoined,stillUnjoined) 
            
    
    def alignDirections(self, superSeg1, superSeg2):    
        overlappingPart=superSeg1.getFirstOverlap(superSeg2)
        superSeg1.makePositive(overlappingPart)
        superSeg2.makePositive(overlappingPart)
        supers=[superSeg1,superSeg2]     
        return supers   
        
    def checkPalindrome(self,segmentList):
        tempInfo=[]
        dePalindromed=[]
        for segment in segmentList:
            tempInfo.append((segment,segment.getOverlap()[0].getName(),segment.getScafStart(), segment.getScafEnd(), segment.getStrand()))
            inOrder=sorted(tempInfo, key=itemgetter(1,2))
        for info in range(len(tempInfo)-1):
            if inOrder[info][2:4]==inOrder[info+1][2:4]:
                if tempInfo[info][4]=='+':
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
        if seg1.getOverlap()[0].getName() == seg2.getOverlap()[0].getName():
            output=[(seg1.getOverlap()[0], 1, seg1.getOverlap()[0].getLength(), seg1.getOverlap()[0].getStrand()),]
            output=self.checkNegatives(output)
            output=SuperSegment(output)
            output.contigs.append(seg1.getContig())
            output.usedScaffolds+=[seg1.getOverlap()[0],]
            return output
                       
        else:
            return self.differentScaffoldBranch(seg1,seg2, joinType)
    
    def diffDirectionBranch(self,seg1,seg2, joinType):
        #print "diffDirectionBranch"
        if seg1.getOverlap()[0].getName() == seg2.getOverlap()[0].getName():
            raise InversionError('''contig %s has a small inversion relative to the reference. please adjust manually.''' % (seg1.getContig().getName()))
        if seg1.getOverlap()[0].getScore() > seg2.getOverlap()[0].getScore():
            newSegTwo = copy.copy(seg2)
            newSegTwo.flipSegment()
            return self.sameDirectionBranch(seg1, newSegTwo, joinType)
        elif (seg1.getOverlap()[0].getScore() == seg2.getOverlap()[0].getScore()) and (seg1.getOverlap()[0].getLength() > seg2.getOverlap()[0].getLength()):
            newSegTwo = copy.copy(seg2)
            newSegTwo.flipSegment()
            return self.sameDirectionBranch(seg1, newSegTwo, joinType)
        else:
            newSegOne = copy.copy(seg1)
            newSegOne.flipSegment()
            return self.sameDirectionBranch(newSegOne, seg2, joinType)
        
        
        
    def differentScaffoldBranch(self,seg1,seg2, joinType):
        #print "differentScaffoldBranch"
        if seg1.getConEnd() < seg2.getConEnd():
            return self.diffScaffCombine(seg1,seg2, joinType) #was notNestedSegsBranch(seg1,seg2)
        else:
            #print "nestedSegs"
            raise NestedSegsError('''Contig %s has nested segments %s and %s''' % (seg1.getContig().getName(), seg1.getName(), seg2.getName()))

    def diffScaffCombine(self,seg1,seg2, joinType):
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
            if joinType == "exclusive":
                output=[(startScaf, 1, startScaf.getLength(), startScaf.getStrand()),]
                output=self.checkNegatives(output)
                output=SuperSegment(output)
                output.contigs.append(contig)
                output.usedScaffolds+=[seg1.getOverlap()[0],seg2.getOverlap()[0]]
            elif joinType == "inclusive":
                output=self.inclusiveSegJoin(seg1,seg2, startScaf.getName())
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
                    
    def inclusiveSegJoin(self,seg1,seg2,envelopingScafName):
        if seg1.getStrand() =='-':
            seg1.flipSegment()
            seg2.flipSegment()
        if seg1.getOverlap()[0].getName() == envelopingScafName:
            enveloper=seg1.getOverlap()[0]
            enveloped=seg2.getOverlap()[0]
        elif seg2.getOverlap()[0].getName() == envelopingScafName:
            enveloper=seg2.getOverlap()[0]
            enveloped=seg1.getOverlap()[0]            
        start=(enveloper,1,seg1.getScafEnd(), enveloper.getStrand())
        middle = (enveloped,seg2.getScafStart(), seg2.getScafEnd(), enveloped.getStrand())
        end = (enveloper,seg1.getScafEnd() + (seg2.getConEnd()-seg1.getConEnd()) ,enveloper.getLength(), seg1.getOverlap()[0].getStrand())
        output=[start,middle,end]
        output=self.checkNegatives(output)
        output=SuperSegment(output)
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
                    

            
           