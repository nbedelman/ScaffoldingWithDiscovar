#!/usr/bin/env python

###take a finished assembly, find areas where parts of scaffolds have been cut
###off, extract those regions, and create new "unplaced" scaffolds with the
###sequences

from Bio import SeqIO
import sys
import csv

lengthFile=sys.argv[1]
assembly=sys.argv[2]

def lengthsDict(tsvFile):
    lengthFile=open(tsvFile, "r")
    outDict={}
    for line in lengthFile:
        scaf=line.split()[0]
        try:
            length=int(line.split()[1])
            outDict[scaf]=length
        except ValueError:
            pass
    return outDict

def findRemovedPieces(lengthTSV,assemblyTSV):
    lengths=lengthsDict(lengthTSV)
    output=[]
    for line in assemblyTSV:
        atts=line.split("\t")
        scaf=atts[2]
        start=int(atts[3])
        end=int(atts[4])
        strand=atts[5]
        if "Hmel" in scaf:
            if start != 1:
                output.append([scaf,1,start-1,strand])
            if end != lengths[scaf]:
                output.append([scaf,end+1,lengths[scaf],strand])
    return output

findRemovedPieces(lengthFile,assembly)
