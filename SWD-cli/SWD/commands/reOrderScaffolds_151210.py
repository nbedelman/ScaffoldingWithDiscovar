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


def runAll(bedDirectory, agpBedFile, originalGenome, discovarAssembly, combineMethod="first", reportDirectory=None):
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
        *fixed.fasta: fasta file with the result of using discovar to join reference scaffolds
    '''
    print "reading contigs"
    rawContigs=readAllContigs(bedDirectory, reportDirectory)
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
        print "writing output"
        print chromosome.getName()
        chromosome.makeAllGroups()
        for group in chromosome.getGroups():
            group.makeSuperScaffolds()
        chromosome.combineGroups(combineMethod)
        chromosome.writeOverviewResults()
        #chromosome.writeFasta(originalGenome, discovarAssembly)
        print "done"
    print "COMPLETED"
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
            base=file.strip(".bed")
            if base not in rejectList:
                con=Contig(rootdir+"/"+file)
                if base in inclusiveList:
                    con.inOrEx = 'i'
                contigs.append(con)
    return contigs
    
def cullSegments(contigList):
    '''small script to loop through each contig and get rid of extraneous alignemtns'''
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
    return combined    

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
        if (not scaffold.getName()=="Ns"):
            try:
                chromDict[scaffold.getChrom()][1].append(scaffold)
            except KeyError:
                chromDict[scaffold.getChrom()]=([],[scaffold])
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
        

chromosomes=runAll("../../../data/fullOverlaps/","../../../data/agpToBed_chroms.bed", \
"../../../data/Hmel2.fa", "../../../data/h_melpomene_clipped_1000.fasta", combineMethod="first", reportDirectory="../../../results_160314/")      


