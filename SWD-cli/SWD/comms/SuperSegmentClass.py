from Bio.Seq import Seq
from Bio import SeqIO
import os
import itertools
import string
import copy
from operator import itemgetter
from ContigClass import *
from ScaffoldClass import *
from PartClass import *
from ErrorClasses import *

class SuperSegment(object):
    def __init__(self, listOfParts):
        self.contigs=[]
        self.scaffolds=[]
        self.partsInOrder=[]
        for part in listOfParts:
            if type(part[0]) == Contig:
                self.contigs.append(part[0])
            elif type(part[0]) == Scaffold:
                self.scaffolds.append(part[0])
            self.partsInOrder.append(Part(part))
        self.envelopers={}
        self.enveloped=[]
        self.usedScaffolds=[]
        self.segIgnores=[]
        self.unusedScaffolds=[]
        self.inclusiveJoin=[]
        self.exclusiveJoin=[]
    def getContigs(self):
        '''return the list of unique contigs used by this superSegment'''
        unique=[]
        uniqueNames=[]
        for contig in self.contigs:
            if contig.getName() not in uniqueNames:
                unique.append(contig)
                uniqueNames.append(contig.getName())
        return unique
    def getScaffolds(self):
        '''return the list of unique scaffolds used by this superSegment'''
        unique=[]
        uniqueNames=[]
        for scaffold in self.scaffolds:
            if scaffold.getName() not in uniqueNames:
                unique.append(scaffold)
                uniqueNames.append(scaffold.getName())
        return unique
    def getPartsInOrder(self):
        return self.partsInOrder
    def getInclusiveJoin(self):
        return self.inclusiveJoin
    def getEnveloped(self):
        return self.enveloped
    def getEnveloperDict(self):
        self.findEnvelopers()
        return self.envelopers
    def getUsedScaffolds(self):
        unique=[]
        uniqueNames=[]
        for scaf in self.usedScaffolds:
            if scaf.getName() not in uniqueNames:
                unique.append(scaf)
                uniqueNames.append(scaf.getName())
        return unique
    def getSegIgnores(self):
        return self.segIgnores
    def getUnusedScaffolds(self):
        unused=[]
        unusedNames=[]
        for scaf in self.unusedScaffolds:
            if scaf.getName() not in unusedNames:
                unused.append(scaf)
                unusedNames.append(scaf.getName())
        return self.unusedScaffolds
    def getTrueIgnores(self):
        ignores=[]
        usedNames=[scaf.getName() for scaf in self.getUsedScaffolds()]
        for ignore in self.getSegIgnores():
            if ignore.getOverlap()[0].getName() not in usedNames:
                ignores.append(ignore)
        return ignores
            
    def reverseOrder(self):
        switched = []
        for item in reversed(self.getPartsInOrder()):
            switched.append(item)
        self.partsInOrder = switched
    
    def flipSuperSegment(self):
        self.reverseOrder()
        for part in self.getPartsInOrder():
            part.flipPart()
            
    def printSuperSeg(self):
        output=''
        for part in self.getPartsInOrder():
            output+=part.printPart()
        return output
    
    def getScore(self):
        scores=[]
        for scaffold in self.getScaffolds():
            scores.append(scaffold.getScore())
        return max(scores)
        
    def getTotalLength(self):
        lengths=[]
        for part in self.getPartsInOrder():
            lengths.append(part.getLength())
        return sum(lengths)
    
    def getFirstCoordinate(self, base):
        sortParts=[]
        for part in range(len(self.getPartsInOrder())):
            if self.getPartsInOrder()[part].getType() == Scaffold:
                sortParts.append((self.getPartsInOrder()[part], self.getPartsInOrder()[part].getScore(), self.getPartsInOrder()[part].getBackboneLength()))
        if base =="first":
            sortParts = sorted(sortParts, key=itemgetter(2), reverse=True)
        elif base=="best":
            sortParts = sorted(sortParts, key=itemgetter(1,2), reverse=True)
        out = sortParts[0][0]
        return out.getBackbone().getStart()
    def isOverlapping(self, otherSuperSegment):
        overlaps=False
        for part in self.getPartsInOrder():
            for otherPart in otherSuperSegment.getPartsInOrder():
                if part.getName() == otherPart.getName():
                    overlaps = True
                    break
        return overlaps
        
    def getFirstOverlap(self, otherSuperSegment, partType=None):
        '''assumes superSegments are oriented in the same direction'''
        thisSuperParts=[]
        otherSuperParts=[]
        for part in self.getPartsInOrder():
            if (not partType) or (part.getType()==partType):
                thisSuperParts.append(part.getName())
        for otherPart in otherSuperSegment.getPartsInOrder():
            if (not partType) or (otherPart.getType()==partType):
                otherSuperParts.append(otherPart.getName())
        for name in thisSuperParts:
            for otherName in otherSuperParts:
                if name== otherName:
                    return name
    
    def getOverlappingIndices(self, overlapName):
        result=[]
        for part in range(len(self.getPartsInOrder())):
            if self.getPartsInOrder()[part].getName() == overlapName:
                result.append(part)
        return result
    
    def lengthBefore(self, backboneName):
        for part in range(len(self.getPartsInOrder())):
            if self.getPartsInOrder()[part].getName() == backboneName:
                index=part
                break #it will stop after first occurrence
        length=0
        for part in range(len(self.getPartsInOrder()[:index])):
            length += self.getPartsInOrder()[part].getLength()
        #length+=self.getPartsInOrder()[index].getStart()
        return length
    
    def lengthAfter(self, backboneName):
        for part in range(len(self.getPartsInOrder())):
            if self.getPartsInOrder()[part].getName() == backboneName:
                index=part #No break: it will continue until the last occurrence
        length=0
        for part in range(index+1,len(self.getPartsInOrder())):
            length += self.getPartsInOrder()[part].getLength()
        #length += (self.getPartsInOrder()[index].getBackboneLength() - self.getPartsInOrder()[index].getEnd())
        return length        
        
    def makePositive(self, backboneName):
        for part in self.getPartsInOrder():
            if part.getName() == backboneName:
                if part.getStrand() == '-':
                    self.flipSuperSegment()
                    return
    
    def alignWithBestScaf(self):
        sortParts=[]
        for part in range(len(self.getPartsInOrder())):
            if self.getPartsInOrder()[part].getType() == Scaffold:
                sortParts.append((part, self.getPartsInOrder()[part].getScore(), self.getPartsInOrder()[part].getBackboneLength()))
        sortParts = sorted(sortParts, key=itemgetter(1,2), reverse=True)
        bestIndex=sortParts[0][0]
        bestPart=self.getPartsInOrder()[bestIndex]
        if bestPart.getStrand() != bestPart.getOriginalStrand():
            self.flipSuperSegment()
        return
    
    def findEnvelopers(self):
        '''figure out if there are scaffolds that lie within other scaffolds
        in this superSegment.'''
        inFinal=[]
        unUsed=[]
        unMatched=[]
        ignored=[ignore.getOverlap()[0].getName() for ignore in self.getTrueIgnores()]
        for part in self.getPartsInOrder():
            if part.getType() == Scaffold:
                inFinal.append(part.getBackbone().getName())
        for scaffold in self.getUsedScaffolds():
            if scaffold.getName() not in inFinal:
                unUsed.append(scaffold.getName())
        for contig in self.getContigs():
            currentConnects=contig.getConnectors()
            missing=[]
            present=[]
            for connect in currentConnects:
                if connect.getName() in unUsed and connect.getName() not in ignored:
                    missing.append(connect.getName())
                elif connect.getName() in inFinal:
                    present.append(connect.getName())
            if missing != []:
                unMatched.append((present,contig.getName(),missing))
                
        matched=self.findMatchingEnvelopes(unMatched, {}, 0)        
        self.enveloped=unUsed
        self.envelopers=matched
         
                
    def findMatchingEnvelopes(self, unMatchedList, matchDict, loopCount):
        unGrounded=[]
        matches=copy.copy(matchDict)
        for item in unMatchedList:
            if item[0] != []:
                ground=item[0][0]
                contig=item[1]
                envelopes=item[2]
                try:
                    matches[ground].append((contig,envelopes+[ground,]))
                    matches[ground][0]+=envelopes
                except KeyError:
                    matches[ground]=[envelopes,(contig,envelopes+[ground,])]
            else:
                unGrounded.append(item)
        stillUnMatched=[]
        for floating in unGrounded:
            contig=floating[1]
            envelopes=floating[2]
            matched=False
            for key in matches.keys():
                contained=matches[key][0]
                for env in envelopes:
                    if env in contained:
                        matched=True
                        for scaf in envelopes:
                            if scaf not in matches[key][0]:
                                matches[key][0].append(scaf)
                        matches[key].append((contig,envelopes))
            if not matched:
                stillUnMatched.append(floating)
        if stillUnMatched != [] and loopCount<10:
            loopCount+=1
            return self.findMatchingEnvelopes(stillUnMatched, matches, loopCount)
        else:
            return matches
                        
        
    
    def hasNewInfo(self, otherSuper):
        unique=False
        otherScaffolds=[]
        for scaf in otherSuper.getUsedScaffolds():
            otherScaffolds.append(scaf.getName())
        for scaffold in self.getScaffolds():
            if scaffold.getName() not in otherScaffolds:
                unique=True
        return unique
        