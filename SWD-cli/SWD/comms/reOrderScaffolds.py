#Take an AGP file, output a bed file that can be used with bedTools to give a fasta file.
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
from ChromosomeClass import *
from SuperSegmentClass import *
from PartClass import *
from GroupClass import *
from ErrorClasses import *

def runAll(bedDirectory, agpBedFile, originalGenome, discovarAssembly, ungroupedChrom='', combineMethod="first", reportDirectory=None, consecutiveOnly=False):
    '''runs the full program.
    Takes:
        bedDirectory: a directory containing a bed file for each discovar contig's alignment to genome
        agpBedFile: a bed file that gives information about where scaffolds are placed on the genome
        originalGenome: a fasta file of the genome
        discovarAssembly: a fasta file of the discovar assembly (a.lines.fasta)
        combineMethod: either "first" or "best". used to place combined scaffolds onto chromosome, based either on
            where the earliest scaffold was originally mapped or where the best (anchored or largest) scaffold was originally mapped
    Returns:
        *coordinates.txt: text files with coordinated of combined scaffolds in form <SCAFFOLD,START,STOP,STRAND>
        *envelopers.txt: text files with information about scaffolds that mapped inside others in form
            ENVELOPING SCAF: <list of enveloped scafs> <list of "evidence" in form (DISCOVAR CONTIG, [SCAFFOLDS COMBINED BY CONTIG])>
        *discoOrder.fasta: fasta file with the result of using discovar to join reference scaffolds
    '''
    print ("reading contigs")
    #If you don't have an "ungrouped" chrom, or you don't want to use it as such, just put in None for ungroupedChrom.
    Contig.ungroupedChrom=ungroupedChrom
    Group.consecutiveOnly=consecutiveOnly
    rawContigs=readAllContigs(bedDirectory, reportDirectory)
    print ("done")
    print ("culling extraneous contig sub-alignments")
    cullSegments(rawContigs)
    print ("done")
    print ("reading scaffolds")
    scaffolds=readScaffold(agpBedFile)
    print ("done")
    print ("finding overlapping scaffolds and contigs")
    for contig in rawContigs:
        contig.findConnectors(scaffolds, "good")
    print ("done")
    print ("combining contig sub-alignments")
    simpContigs=combineSegments(rawContigs)
    print ("done")
    print ("finding connections after simplifying")
    for contig in simpContigs:
        contig.findConnectors(scaffolds, "combined")
    print ("grouping by Chromosome")
    chromDict=groupPiecesByChromosome(simpContigs,scaffolds)
    chromosomes=makeChromosomes(chromDict)
    print ("done")
    print ("grouping within Chromosome and making super scaffolds")
    usedUngroupedScafs=[]
    unGrouped=False
    for chromosome in chromosomes:
        if chromosome.getName()!=Contig.ungroupedChrom:
            print ("writing output")
            print (chromosome.getName())
            chromosome.makeAllGroups()
            print ("made all groups")
            #changed this on 6/14/17 so I can recycle parts of groups that didn't really belong.
            #numGroups=len(chromosome.getGroups())
            for group in chromosome.getGroups():#range(numGroups):
                group.makeSuperScaffolds() #chromosome.getGroups()[group].makeSuperScaffolds()
            chromosome.combineGroups(combineMethod)
            usedUngroupedScafs+=chromosome.getUsedUngroupedScafs()
            print("combined groups")
            chromosome.writeCoordinates(chromosome.getSuperScaffolds(),overwritten=False)
            chromosome.writeCoordinates(chromosome.getOverwrittenSupers(), overwritten=True)
            chromosome.writeOverviewResults()
            print("wrote results")
            chromosome.writeFasta(originalGenome, discovarAssembly)
            print ("done")
        else:
            unGrouped=chromosome
    if Contig.ungroupedChrom and unGrouped:
        unGrouped.removeScaffolds(usedUngroupedScafs)
        print ("writing output")
        print (unGrouped.getName())
        unGrouped.makeAllGroups()
        print ("made all groups")
        numGroups=len(unGrouped.getGroups())
        for group in range(numGroups):
            unGrouped.getGroups()[group].makeSuperScaffolds()
        unGrouped.combineGroups(combineMethod)
        print("combined groups")
        unGrouped.writeCoordinates(chromosome.getSuperScaffolds(),overwritten=False)
        unGrouped.writeCoordinates(chromosome.getOverwrittenSupers(), overwritten=True)
        unGrouped.writeOverviewResults()
        print("wrote results")
        unGrouped.writeFasta(originalGenome, discovarAssembly)
    print ("done")
    print ("COMPLETED")
    return chromosomes


def readAllContigs(directory, reportDirectory):
    '''takes a directory with bed files of contig alignments
    returns a list of Contig objects'''
    rejectList=[]
    inclusiveList=[]
    if reportDirectory:
        for subdir, dirs, files in os.walk(reportDirectory):
            for file in files:
                if "report" in file:
                    data=read_csv(reportDirectory+file)
                    for line in data:
                        if line[6].lower()=="x":
                            rejectList.append(line[0])
                            write_csv("reject_report.csv",[line,])
                        elif line[6].lower()=="i":
                            inclusiveList.append(line[0])
    contigs=[]
    rootdir = directory
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if ".bed" in file:
                base=file.strip(".bed")
                if base not in rejectList:
                    con=Contig(rootdir+"/"+file)
                    if base in inclusiveList:
                        con.inOrEx = 'i'
                    contigs.append(con)
    return contigs

def cullSegments(contigList):
    '''small script to loop through each contig and get rid of extraneous aligments'''
    for i in contigList:
        i.cullSegments()
    return

def combineSegments(contigList):
    '''small script to loop through each contig and combine alignments that are very close to one another'''
    combined=[]
    for i in contigList:
        try:
            i.combineSegments()
            combined.append(i)
        except NotInformativeError:
            pass
        except NoSegmentError:
            pass
    return combined

def readScaffold(bedFile):
    '''take a bed file
    return a list of objects of class scaffold, in the order they are in the file'''
    #6/27/27 using a list format is making the findOverlaps step extremely slow.
    #I think a dictionary would make it faster.
    #bedList=[]
    #with open(bedFile,"r") as f:
    #    for line in f:
    #        scaf=Scaffold(line)
    #        if not scaf.getName()=='Ns':
    #            bedList.append(scaf)
    #return bedList
    bedDict={}
    with open(bedFile, "r") as f:
        for line in f:
            scaf=Scaffold(line)
            if 'Ns' not in scaf.getName():
                try:
                    bedDict[scaf.getChrom()][0].append(scaf)
                except KeyError:
                    bedDict[scaf.getChrom()]=[[scaf,],[]]
            else:
                try:
                    bedDict[scaf.getChrom()][1].append(scaf)
                except KeyError:
                    bedDict[scaf.getChrom()]=[[],[scaf,]]
    return bedDict


def groupPiecesByChromosome(contigs, scaffolds):
    chromDict={}
    for contig in contigs:
        try:
            chromDict[contig.getChrom()][0].append(contig)
        except KeyError:
            chromDict[contig.getChrom()]=[[contig,],[]]
    for scaffold in scaffolds.keys():
        #want N's as well as non-Ns here.
        try:
            chromDict[scaffold][1]=scaffolds[scaffold][0]+scaffolds[scaffold][1]
        except KeyError:
            chromDict[scaffold]=[[],scaffolds[scaffold][0]+scaffolds[scaffold][1]]
        #updated 6/27/17 because scaffolds is now a dict.
        #if (not scaffold.getName()=="Ns"):
        #    try:
        #        chromDict[scaffold.getChrom()][1].append(scaffold)
        #    except KeyError:
        #        chromDict[scaffold.getChrom()]=([],[scaffold])
    return chromDict


def makeChromosomes(chromDict):
    chromosomes=[]
    for key in chromDict.keys():
        chromosomes.append(Chromosome(key, chromDict[key][0], chromDict[key][1]))
    return chromosomes

def read_csv(file, header=None):
    '''
    Takes a csv file,
    outputs a list, where each element is a list which contains the cells of a single row as strings.
    '''
    data = []
    reader = csv.reader(open(file, 'rU'))
    for row in reader:
        data.append(row)
    if header == True:
        return data[1:]
    else:
        return data

def write_csv(file, asequence, header=None):
    fp =  open(file, 'a+')
    a = csv.writer(fp, delimiter=',')
    if header:
        a.writerows(header)
    a.writerows(asequence)
    fp.close()
