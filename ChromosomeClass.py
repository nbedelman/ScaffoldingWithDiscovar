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
        return Group(finishedLists[0], finishedLists[1], self)
        
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
        totalUsed=[]
        totalUnused=[]
        for group in self.getGroups():
            if group.getFullGroupJoin() != []:
                allSuperScafs.append(group.getFullGroupJoin())
                totalUsed+=group.getFullGroupJoin().getEnveloped()
                totalUsed+=[scaf.getName() for scaf in group.getFullGroupJoin().getUsedScaffolds()]
                if group.getFullGroupJoin().getTrueIgnores() != [] or group.getFullGroupJoin().getUnusedScaffolds() != []:
                    unGroupedNames=[scaf.getName() for scaf in self.getUngroupedScaffolds()]
                    allUnused=[ignore.getOverlap()[0] for ignore in group.getFullGroupJoin().getTrueIgnores()] + \
                                group.getFullGroupJoin().getUnusedScaffolds()
                    for ignore in allUnused:
                        if (ignore.getName() not in unGroupedNames) and (not self.checkIfUsed(self.unGroupedScaffolds, ignore)):
                            totalUnused.append(ignore)
                            unGroupedNames.append(ignore.getName())
        
            else:
                for scaffold in group.getScaffoldList():
                    if (not self.checkIfUsed(self.unGroupedScaffolds, scaffold)) and (scaffold.getName() not in totalUsed):
                        totalUnused.append(scaffold)
        
        for unUsed in totalUnused:
            if not unUsed.getName() in totalUsed:
                self.unGroupedScaffolds.append(unUsed)
                
        for scaffold in self.getUngroupedScaffolds():
            superScaf=SuperSegment([(scaffold,1,scaffold.getLength(),scaffold.getStrand()),])
            allSuperScafs.append(superScaf)
        sortable=[]
        for superScaf in allSuperScafs:
            sortable.append((superScaf, superScaf.getFirstCoordinate(base)))
        sortable=sorted(sortable, key=itemgetter(1))
        onlyScafs=[sup[0] for sup in sortable]
        self.superScaffolds=onlyScafs
        
       
    def writeOverviewResults(self):
        f=open(self.getName()+"coordinates.txt", "w")
        env=open(self.getName()+"envelopers.txt","w")
        groupCount=0
        coords=[]
        envs=[]
        for superScaf in self.getSuperScaffolds():
            if superScaf.printSuperSeg() not in coords:
                f.write("GROUP NUMBER: " + str(groupCount) + "\n")
                f.write(superScaf.printSuperSeg())
                groupCount+=1
                coords.append(superScaf.printSuperSeg())
        for group in self.getGroups():
            fullJoin=group.getFullGroupJoin()
            if fullJoin !=[] and group.printGroupEnvelopers() not in envs:
                env.write(group.printGroupEnvelopers())
                envs.append(group.printGroupEnvelopers())
        f.close
        env.close
        
    def writeFasta(self, chromOrScaf, originalScaffolds, originalDiscoContigs):
        originalScafs=SeqIO.index(originalScaffolds, "fasta")
        originalDiscos=SeqIO.index(originalDiscoContigs, "fasta")
        groupCount=0
        coords=[]
        newScafs=[]
        for superScaf in self.getSuperScaffolds():
            newRecord=SeqIO.SeqRecord(seq='')
            theseCoords=superScaf.printSuperSeg()
            records=theseCoords.split("\n")[:-1]
            usedParts=[]
            for record in records:
                items=record.split(",")
                name=items[0]
                start=int(items[1])-1
                end=int(items[2])-1
                direction=items[3]
                try:
                    original=originalScafs[name]
                except KeyError:
                    original=originalDiscos[name]
                if direction == '-':
                    original=original.reverse_complement()
                newPart=original[start:end]
                newRecord=newRecord+newPart
                usedParts.append(name)
            newRecord.id=self.getName()+"NewScaf"+str(groupCount)
            newRecord.description=str(usedParts)
            newScafs.append(newRecord)
            groupCount+=1
        SeqIO.write(newScafs,open(self.getName()+"fixed.fasta", "w"), "fasta")
            
            
