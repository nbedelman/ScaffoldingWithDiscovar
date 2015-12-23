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
from GroupClass import *
from ErrorClasses import *

class Chromosome(object):
    def __init__(self, chromName, contigList, scafList):
        self.name=chromName
        self.scaffoldList=scafList
        self.contigList=contigList
        self.groups=[]
        self.unGroupedScaffolds=[]
        self.superScaffolds=[]
    def getName(self):
        return self.name
    def getScaffoldList(self):
        return self.scaffoldList
    def getContigList(self):
        return self.contigList
    def getGroups(self):
        return self.groups
    def getUngroupedScaffolds(self):
        return self.unGroupedScaffolds 
    def getSuperScaffolds(self):
        return self.superScaffolds       
    def makeAllGroups(self):
        allContigs=copy.copy(self.contigList)
        usedContigs=[]
        usedScaffolds=[]
        for contig in allContigs:
            if not self.checkIfUsed(usedContigs, contig):
                group=self.makeGroup(contig)
                if group.getScaffoldList() != []:
                    usedContigs += group.getContigList()
                    usedScaffolds+=group.getScaffoldList()
                    self.groups.append(group)
        for scaffold in self.getScaffoldList():
            if not self.checkIfUsed(usedScaffolds, scaffold):
                self.unGroupedScaffolds.append(scaffold)
    
    def checkIfUsed(self,aList, anObject):
        used=False
        for item in aList:
            if item.getName() == anObject.getName():
                used=True
                break
        return used    
    
    def makeGroup(self,contig):
        contig_list=[contig,]
        scaffold_list = copy.copy(contig.getConnectors())
        listNow=[scaf.getName() for scaf in scaffold_list]
        finishedLists=(self.backAndForth(contig_list,scaffold_list))
        return Group(finishedLists[0], finishedLists[1])
        
    def backAndForth(self,contigList,scaffoldList):
        startContigs=len(contigList)
        startScaffolds=len(scaffoldList)
        newContigs=[]
        newScaffolds=[]
        for scaffold in scaffoldList:
            for overlappingContig in scaffold.getOverlaps():
                if not self.checkIfUsed(newContigs, overlappingContig):
                    newContigs.append(overlappingContig)
        for contig in contigList:
            for overlappingScaffold in contig.getConnectors():
                if not self.checkIfUsed(newScaffolds, overlappingScaffold):   
                    newScaffolds.append(overlappingScaffold)       
        if startContigs == len(newContigs) and startScaffolds==len(newScaffolds):
            return (newContigs,newScaffolds)
        else:
            return self.backAndForth(newContigs,newScaffolds)
            
    def combineGroups(self, base):
        #parallelize form, and put all the segments in a single place
        allSuperScafs=[]
        for group in self.getGroups():
            if group.getFullGroupJoin() != []:
                allSuperScafs.append(group.getFullGroupJoin())
                if group.getFullGroupJoin().getTrueIgnores() != []:
                    unGroupedNames=[scaf.getName() for scaf in self.getUngroupedScaffolds()]
                    for ignore in group.getFullGroupJoin().getTrueIgnores():
                        if ignore.getOverlap()[0].getName() not in unGroupedNames:
                            self.unGroupedScaffolds.append(ignore.getOverlap()[0])
                            unGroupedNames.append(ignore.getOverlap()[0].getName())
            else:
                for scaffold in group.getScaffoldList():
                    self.unGroupedScaffolds.append(scaffold)
        for scaffold in self.getUngroupedScaffolds():
            superScaf=SuperSegment([(scaffold,1,scaffold.getLength(),scaffold.getStrand()),])
            allSuperScafs.append(superScaf)
        sortable=[]
        for superScaf in allSuperScafs:
            sortable.append((superScaf, superScaf.getFirstCoordinate(base)))
        sortable=sorted(sortable, key=itemgetter(1))
        onlyScafs=[sup[0] for sup in sortable]
        self.superScaffolds=onlyScafs
        
       
    def writeResults(self):
        f=open(self.getName()+"coordinates.txt", "w")
        env=open(self.getName()+"envelopers.txt","w")
        groupCount=0
        for superScaf in self.getSuperScaffolds():
            f.write("GROUP NUMBER: " + str(groupCount) + "\n")
            f.write(superScaf.printSuperSeg())
            groupCount+=1
        for group in self.getGroups():
            fullJoin=group.getFullGroupJoin()
            if fullJoin:
                env.write(group.printGroupEnvelopers())
        f.close
        env.close
