import csv
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


def write_csv(file, asequence, header=None):
    fp =  open(file, 'w')
    a = csv.writer(fp, delimiter=',')
    if header:
        a.writerows(header)
    a.writerows(asequence)
    fp.close()

class Chromosome(object):
    def __init__(self, chromName, contigList, scafList):
        self.name=chromName
        self.scaffoldList=scafList
        self.contigList=contigList
        self.groups=[]
        self.unGroupedScaffolds=[]
        self.superScaffolds=[]
        self.usedUngroupedScafs=[]
        self.overwrittenSupers=[]
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
    def getUsedUngroupedScafs(self):
        return self.usedUngroupedScafs
    def getOverwrittenSupers(self):
        return self.overwrittenSupers
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

    def removeScaffolds(self, listOfScafNames):
        currentScafs=self.getScaffoldList()
        self.scaffoldList=[]
        for s in currentScafs:
            if s.getName() not in listOfScafNames:
                self.scaffoldList.append(s)

    def checkIfUsed(self,aList, anObject):
        used=False
        for item in aList:
            if item.getName() == anObject.getName():
                used=True
                break
        return used

    def makeGroup(self,contig):
        if self.getName() != Contig.ungroupedChrom:
            contig_list=[contig,]
            scaffold_list = copy.copy(contig.getConnectors())
            finishedLists=(self.backAndForth(contig_list,scaffold_list))
            return Group(finishedLists[0], finishedLists[1], self)
        #we don't want to connect placed scaffolds to unplaced ones in the unplaced scaffold chrom. The only thing we want to do is extend if possible.
        else:
            if len(contig.getConnectors())==1:
                return Group([contig,],contig.getConnectors(),self)
            else:
                return Group([contig,],[],self)

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
                #Need to make sure we're not connecting across different chromosomes. Don't add scaffolds outside the current chrom/ungroupedChroms
                if (overlappingScaffold.getChrom() == self.getName()) or (overlappingScaffold.getChrom() == Contig.ungroupedChrom):
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
            if scaffold.getChrom() == self.getName():
                superScaf=SuperSegment([(scaffold,0,scaffold.getLength(),scaffold.getStrand()),])
                allSuperScafs.append(superScaf)

        self.arrangeSuperScafs(allSuperScafs, base)

        #Record if any of the ungroupedChrom scaffolds were included. Only use this if actually using an ungroupedChrom!
        if Contig.ungroupedChrom:
            for scaffold in totalUsed:
                if Contig.ungroupedChrom in scaffold:
                    self.usedUngroupedScafs.append(scaffold)

    def arrangeSuperScafs(self,superScafList, base):
        sortable=[]
        include=[]
        #first, sort the superScaffolds by whatever metric is specified by base
        for superScaf in superScafList:
            firstCoord=superScaf.getSortCoordinates(base)[0]
            lastCoord=superScaf.getSortCoordinates(base)[1]
            sortable.append((superScaf,firstCoord, lastCoord))
        #want to sort by start position, earliest first, and then size, largest first. The following two lines accomplish this.
        sortable=sorted(sortable, key=itemgetter(2), reverse=True)
        sortable=sorted(sortable, key=itemgetter(1))

        #Then, get rid of scaffolds that have been overwritten and put them in a separate variable. Some of these are going to be duplicates; I don't know why.
        include.append(sortable[0])
        for sortScaf in sortable[1:]:
            lastIn=include[-1]
            if (lastIn[1]<sortScaf[1]) and (lastIn[2]>sortScaf[2]):
                self.overwrittenSupers.append(sortScaf[0])
            elif (lastIn[2]==sortScaf[2]) or (lastIn[1]==sortScaf[1]):
                pass
            elif lastIn[2]>sortScaf[1]:
                if (lastIn[0].getJoinType()=="ext") and ("Ns" in sortScaf[0].getPartsInOrder()[0].getName()):
                    sortScaf[0].getPartsInOrder()[0].start+=lastIn[2]-sortScaf[1]+1
                    include.append(sortScaf)
                elif (sortScaf[0].getJoinType()=="ext") and ("Ns" in lastIn[0].getPartsInOrder()[0].getName()):
                    lastIn[0].getPartsInOrder()[0].end-=lastIn[2]-sortScaf[1]-1
                    include.append(sortScaf)
                else:
                    include.append(sortScaf)
            else:
                include.append(sortScaf)
        onlyScafs=[sup[0] for sup in include]
        self.superScaffolds=onlyScafs



    def writeCoordinates(self, superScafList, overwritten):
        if overwritten:
            nameBase=self.getName()+"_overwrittenScaffolds"
            f=open(nameBase+".txt", "w")
        else:
            nameBase=self.getName()+"_"
            f=open(nameBase+"coordinates.txt", "w")
        groupCount=0
        coords=[]
        for superScaf in superScafList:
            if superScaf.printSuperSeg() not in coords:
                f.write(nameBase+"number_" + str(groupCount) + "\n")
                f.write(superScaf.printSuperSeg()+"\n")
                groupCount+=1
                coords.append(superScaf.printSuperSeg())
        f.close()

    def writeOverviewResults(self):
        env=open(self.getName()+"_envelopers.txt","w")
        envs=[]
        reps=[]
        for group in self.getGroups():
            fullJoin=group.getFullGroupJoin()
            if fullJoin !=[]:# and group.printGroupEnvelopers() not in envs:
                env.write(group.printGroupEnvelopers())
                envs.append(group.printGroupEnvelopers())
                reps+=(group.getContigReports())

        env.close()
        write_csv(self.getName()+"report.csv",reps)

    def writeFasta(self, originalScaffolds, originalDiscoContigs):
        originalScafs=SeqIO.index(originalScaffolds, "fasta")
        originalDiscos=SeqIO.index(originalDiscoContigs, "fasta")
        groupCount=0
        coords=[]
        newScafs=[]
        for superScaf in self.getSuperScaffolds():
            if superScaf.printSuperSeg() not in coords:
                newRecord=SeqIO.SeqRecord(seq='')
                theseCoords=superScaf.printSuperSeg()
                records=theseCoords.split("\n")[:-1]
                usedParts=[]
                for record in records:
                    items=record.split(",")
                    name=items[0]
                    if int(items[1])<0:
                        start=0
                    else:
                        start=int(items[1])
                    end=int(items[2])
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
            coords.append(superScaf.printSuperSeg())
        SeqIO.write(newScafs,open(self.getName()+"discoOrder.fasta", "w"), "fasta")
