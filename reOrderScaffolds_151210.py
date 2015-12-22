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
# for chromosome in chromosomes:
#     if chromosome.getName() == "Hmel204":
#         Hmel204 = chromosome
# for item in Hmel203.getGroups():
#     scafList=[scaf.getName() for scaf in item.getScaffoldList()]
#     if "Hmel204013" in scafList:
#         group=item  

#superScaffolds=[]
#contig=group.getContigList()[0]
#
#segments=[]
#smallIgnores=[]
#for seg in contig.getCombinedSegments():
#    if seg.getLength() < 1000:
#        smallIgnores.append(copy.copy(seg))
#    else:
#        segments.append(copy.copy(seg))
#orderedSegs=group.orderSegs(segments)
#palindromeChecked=group.checkPalindrome(orderedSegs)
#
#
#segmentList=palindromeChecked
#joinedSegments=[]
#ignoredSegs=[]
#
#joinedSegments=group.joinEachSegment(palindromeChecked, [],smallIgnores)
#fullContigJoin=group.joinSuperScaffolds(joinedSegments,[])
#if fullContigJoin != []:
#    if fullContigJoin.getContigs() == []:
#        fullContigJoin.contigs.append(contig)
#    superScaffolds.append(fullContigJoin)
        

# superScaffolds=[]
# for contig in group.getContigList():
#     segments=[]
#     smallIgnores=[]
#     for seg in contig.getCombinedSegments():
#         if seg.getLength() < 1000:
#             smallIgnores.append(copy.copy(seg))
#         else:
#             segments.append(copy.copy(seg))
#     orderedSegs=group.orderSegs(segments)
#     palindromeChecked=group.checkPalindrome(orderedSegs)
#     joinedSegments=group.joinEachSegment(palindromeChecked, [],smallIgnores)
#     fullContigJoin=group.joinSuperScaffolds(joinedSegments,[])
#     if fullContigJoin != []:
#         if fullContigJoin.getContigs() == []:
#             fullContigJoin.contigs.append(contig)
#         superScaffolds.append(fullContigJoin) 
# fullGroupJoin=group.joinSuperScaffolds(superScaffolds,[])
# if fullGroupJoin != []:
#     fullGroupJoin.alignWithBestScaf()
   
               
                      
                             
                                           
                
#joinedSupers=[]                       
#notJoined=[]                                        
#index=0
#toJoin=superScaffolds[index]        
#notJoined.append(superScaffolds[index])           
#              
        
#        
#fullGroupJoin=self.joinSuperScaffolds(superScaffolds,[])
#if fullGroupJoin != []:
#    fullGroupJoin.alignWithBestScaf()
#    fullGroupJoin.segIgnores=trueIgnores
#self.fullGroupJoin=fullGroupJoin
