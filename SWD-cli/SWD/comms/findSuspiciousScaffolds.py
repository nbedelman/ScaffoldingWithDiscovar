#!/bin/env python

import csv
import sys

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

def compareLengths(tsvFile, lengthDict):
    tsv=open(tsvFile, "r")
    outDict={}
    for line in tsv:
        try:
            atts=line.split()
            scaf=atts[2]
            start=int(atts[3])
            stop=int(atts[4])
            length=(stop-start)-1
            try:
                diff=lengthDict[scaf]-length
                outDict[scaf]=[lengthDict[scaf],length,diff,float(diff)/lengthDict[scaf]]
            except KeyError:
                pass
        except ValueError:
            pass
    return outDict


def extractProblemScaffs(comparisonDict):
    problemScaffs={}
    for scaf in comparisonDict:
        diff=comparisonDict[scaf][2]
        pct=comparisonDict[scaf][3]
        if (diff>10000) and (pct > .1):
            problemScaffs[scaf]=comparisonDict[scaf]
    return problemScaffs

def write_csv(outFile, comparisonDict,):
    fp =  open(outFile, 'w')
    a = csv.writer(fp, delimiter=',')
    a.writerows([["Scaffold","originalLength","newLength","Difference","PCT"]])
    comparisonSeq=[]
    for entry in comparisonDict:
        entry_as_seq=[entry,]+comparisonDict[entry]
        comparisonSeq.append(entry_as_seq)
    a.writerows(comparisonSeq)
    fp.close()

#Run the functions#

lengthFile=sys.argv[1]
coordinateFile=sys.argv[2]

lengths=lengthsDict(lengthFile)
comparisonDict=compareLengths(coordinateFile,lengths)
problemScaffs=extractProblemScaffs(comparisonDict)
write_csv("lengthComparison.csv",problemScaffs)
