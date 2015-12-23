#Take an AGP file, output a bed file that can be used with bedTools to give a fasta file.
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
from ChromosomeClass import *
from SuperSegmentClass import *
from PartClass import *
from GroupClass import *
from ErrorClasses import *


def runAll(bedDirectory, agpBedFile):
    print "reading contigs"
    rawContigs=readAllContigs(bedDirectory)
    print "done"
    print "culling extraneous contig sub-alignments"
    cullSegments(rawContigs)
    print "done"
    print "reading scaffolds"
    scaffolds=readScaffold(agpBedFile)
    print "done"
    print "finding overlapping scaffolds and contigs"
    for contig in rawContigs:
        contig.findConnectors(scaffolds)
    print "done"
    print "combining contig sub-alignments"
    simpContigs=combineSegments(rawContigs)
    print "done"
    print "grouping by Chromosome"
    chromDict=groupPiecesByChromosome(simpContigs,scaffolds)
    chromosomes=makeChromosomes(chromDict)
    print "done"
    print "grouping within Chromosome and making super scaffolds"
    for chromosome in chromosomes:
        print chromosome.getName()
        chromosome.makeAllGroups()
        for group in chromosome.getGroups():
            group.makeSuperScaffolds()
    print "done"
    return chromosomes
    

def readAllContigs(directory):
    contigs=[]
    rootdir = directory
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            contigs.append(Contig(rootdir+"/"+file))
    return contigs
    
def cullSegments(contigList):
    for i in contigList:
        i.cullSegments()
    return
        
def combineSegments(contigList):
    combined=[]
    for i in contigList:
        try: 
            i.combineSegments()
            combined.append(i)
        except NotInformativeError:
            pass
    return combined    

def intersectsTocontigList(intersectFile):
    out='cp '
    contigs=[]
    with open(intersectFile, "r") as f:
        for line in f:
            atts=line.split("\t")
            if int(atts[-1])>1:
                contigs.append(atts[3])
    for i in contigs:
        out+=i+".bed "
    out += "../../overlaps"
    with open("../../copyCMD.sh","w") as o:
        o.write(out)
    return contigs
                

def readScaffold(bedFile):
    '''take a bed file
    return a list of objects of class seqType (scaffolds or contigs), in the order they are in the file'''
    bedList=[]
    with open(bedFile,"r") as f:
        for line in f:
            scaf=Scaffold(line)
            if not scaf.getName()=='Ns':
                bedList.append(scaf)
    return bedList
       
def groupPiecesByChromosome(contigs, scaffolds):
    chromDict={}
    for contig in contigs:
        try:
            chromDict[contig.getChrom()][0].append(contig)
        except KeyError:
            chromDict[contig.getChrom()]=([contig,],[])
    for scaffold in scaffolds:
        if (not scaffold.getChrom()=='Hmel200') and (not scaffold.getName()=="Ns"):
            try:
                chromDict[scaffold.getChrom()][1].append(scaffold)
            except KeyError:
                pass
    return chromDict

def makeChromosomes(chromDict):
    chromosomes=[]
    for key in chromDict.keys():
        chromosomes.append(Chromosome(key, chromDict[key][0], chromDict[key][1]))
    return chromosomes
        
chromosomes=runAll("simulatedData/fullOverlaps/","simulatedData/agpToBed_chroms.bed")      

for chromosome in chromosomes:
    print chromosome.getName()
    chromosome.combineGroups("first")
    chromosome.writeResults()
#
for chromosome in chromosomes:
    if chromosome.getName() == "Hmel201":
         Hmel201 = chromosome
for item in Hmel201.getGroups():
    scafList=[scaf.getName() for scaf in item.getScaffoldList()]
    if "Hmel201019" in scafList:
         group=item
#
#if joinedSupers == [] and len(superScaffolds) > 0:
#            joinedSupers=superScaffolds[0]
#            return self.joinSuperScaffolds(superScaffolds[1:],joinedSupers)
#elif len(superScaffolds) == 0:
#    return joinedSupers
#else:
#    notJoined=[]
#    for index in range(len(superScaffolds)):
#        toJoin=superScaffolds[index]
#        if toJoin.isOverlapping(joinedSupers):
#            if toJoin.hasNewInfo(joinedSupers):
#                joinedSuperScafs=joinedSupers.getUsedScaffolds()
#                joinedSuperContigs=joinedSupers.getContigs()
#                
#                overlappingPart=toJoin.getFirstOverlap(joinedSupers)
#                toJoin.makePositive(overlappingPart)
#                joinedSupers.makePositive(overlappingPart)
#                supers=[joinedSupers,toJoin]
#                
#                joinedSupersStart=joinedSupers.lengthBefore(overlappingPart)
#                toJoinStart=toJoin.lengthBefore(overlappingPart)
#                distFromStart=[joinedSupersStart,toJoinStart]
#                
#                joinedSupersEnd=joinedSupers.lengthAfter(overlappingPart)
#                toJoinEnd=toJoin.lengthAfter(overlappingPart)                        
#                distFromEnd=[joinedSupersEnd,toJoinEnd]
#                
#                furthestFromStart=distFromStart.index(max(distFromStart))
#                startSuper=supers[furthestFromStart]
#                
#                furthestFromEnd=distFromEnd.index(max(distFromEnd))
#                endSuper=supers[furthestFromEnd]
#                
#                startOverlappingIndex=startSuper.getOverlappingIndices(overlappingPart)[0]
#                endOverlappingIndex=endSuper.getOverlappingIndices(overlappingPart)[0]
#                
#                newStart=[]
#                for part in range(len(startSuper.getPartsInOrder()[:startOverlappingIndex])):
#                    newStart.append(startSuper.getPartsInOrder()[part].exportPart())
#                    
#                newEnd=[]
#                try:
#                    for part in range(endOverlappingIndex+1, len(endSuper.getPartsInOrder())):
#                        newEnd.append(endSuper.getPartsInOrder()[part].exportPart())
#                except IndexError:
#                    newEnd=[]
#                
#                
#                startOverlapPart=startSuper.getPartsInOrder()[startOverlappingIndex]
#                endOverlapPart=endSuper.getPartsInOrder()[endOverlappingIndex]
#                overlapBackbone=startOverlapPart.getBackbone()
#                
#                newOverlapPart=[(overlapBackbone,max(startOverlapPart.getStart(), endOverlapPart.getStart()), min(startOverlapPart.getEnd(), endOverlapPart.getEnd()),'+'),]
#                joinedSupers=newStart+newOverlapPart+newEnd
#                joinedSupers=self.checkNegatives(joinedSupers)
#                joinedSupers = SuperSegment(joinedSupers)
#                if toJoin.getContigs() != []:
#                    joinedSupers.contigs.append(toJoin.contigs[0])
#                joinedSupers.contigs+=joinedSuperContigs
#                joinedSupers.usedScaffolds+=toJoin.getUsedScaffolds()
#                joinedSupers.usedScaffolds+=joinedSuperScafs
#                
#                try:
#                    unjoined=notJoined+superScaffolds[index+1:]
#                except IndexError:
#                    unjoined=notJoined
#                return self.joinSuperScaffolds(unjoined, joinedSupers)
#            else:
#                if toJoin.getContigs() != []:
#                    joinedSupers.contigs.append(toJoin.contigs[0])
#                joinedSupers.usedScaffolds+=toJoin.getUsedScaffolds()
#                try:
#                    unjoined=notJoined+superScaffolds[index+1:]
#                except IndexError:
#                    unjoined=notJoined
#                return self.joinSuperScaffolds(unjoined, joinedSupers) 
#        else:
#            notJoined.append(superScaffolds[index]) 
##If there are superScaffolds that could not be joined, it's because we have something like
##scaffolds 1 and 2 both were placed inside 3, and the superScaffold that links 1 and 2 together
##therefore doesn't overlap with the big combined one
#for supScaf in notJoined:
#    if supScaf.getContigs() != []:
#        joinedSupers.contigs.append(supScaf.getContigs()[0])
#    joinIgnores=[seg.getOverlap()[0].getName() for seg in joinedSupers.getTrueIgnores()]
#    for used in supScaf.getUsedScaffolds():
#        if used.getName() not in joinIgnores:
#            joinedSupers.usedScaffolds.append(used)
#return joinedSupers   