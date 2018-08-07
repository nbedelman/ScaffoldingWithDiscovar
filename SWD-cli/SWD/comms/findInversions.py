#!/usr/bin/env python

#findInversions.py
#look for inversions between a DISCOVAR genome and a reference genome. Only consider those deletions found WITHIN a single reference scaffold.
import sys
import csv
from Bio.Seq import Seq
from Bio import SeqIO
import os
import itertools
import string
import copy
from operator import itemgetter


# sys.path.insert(0, '/Users/nbedelman/Documents/Mallet_Lab/referenceScaffolding/ScaffoldingWithDiscovar/SWD-cli/SWD/comms/')


from ScaffoldClass import *
from ContigClass import *
from SegmentClass import *
from ChromosomeClass import *
from SuperSegmentClass import *
from PartClass import *
from InversionCandidateClass import *
from ErrorClasses import *
from reOrderScaffolds import read_csv, write_csv, groupPiecesByChromosome, makeChromosomes, readScaffold, combineSegments, cullSegments, readAllContigs


#agpBedFile="/Users/nbedelman/Documents/Mallet_Lab/referenceScaffolding/Hmel3/Hmel2_ordered_simpleName.fa.bed"
#agpBedFile="/Users/nbedelman/Documents/Mallet_Lab/18Genomes/Genomic_Analysis/correctScaffs_fixed.bed"
# agpBedFile="/Users/nbedelman/Documents/Mallet_Lab/referenceScaffolding/refGenomeChrom.fasta.bed"
# bedDirectory="/Users/nbedelman/Documents/Mallet_Lab/referenceScaffolding/inversionBeds"
# bpFile="/Users/nbedelman/Desktop/inversionCandidates/breakPoints.bed"

agpBedFile=sys.argv[1]
bedDirectory=sys.argv[2]
outDir=sys.argv[3]
pctThreshold=float(sys.argv[4])
minLength=int(sys.argv[5])

scaffolds=readScaffold(agpBedFile)
rawContigs=readAllContigs(bedDirectory,None)
for contig in rawContigs:
    contig.findConnectors(scaffolds, "good")
cullSegments(rawContigs)
simpContigs=combineSegments(rawContigs, multiScafs=False)
for contig in simpContigs:
    contig.findConnectors(scaffolds, 'combined')
o=open(outDir+"/allBreakPoints.bed","w")
for i in simpContigs:
    inv=InversionCandidate(i)
    inv.findInversion(pctThreshold,minLength)
    if inv.isInversion:
        inv.outputBed(inv.getSegments(),'''%s/%s.bed''' % (outDir,inv.getContig().getName()))
        o.write(inv.getBreakPoints())
o.close()
