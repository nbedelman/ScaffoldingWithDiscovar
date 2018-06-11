#!/bin/env python

#coordinatesToTSV.py
#Some aspects of this script are specific to Dryas!
import sys

coordinateFile=sys.argv[1]
outFile=sys.argv[2]

def coordinatesToTSV(coordinateFile,outFile):
    c=open(coordinateFile,"r")
    o=open(outFile, "w")
    header='''%s\t%s\t%s\t%s\t%s\t%s\n''' % ("Chrom","Group","Scaffold","Start","Stop","Strand")
    o.write(header)
    currentChrom=''
    currentGroup=''
    for line in c:
        if "number" in line:
            name=line.split("_")
            #this one works for dryas
            # currentChrom=name[0]+"_"+name[1]
            # currentGroup=name[3].strip("\n")
            #This one works for Hmel
            currentChrom=name[0]
            currentGroup=name[2].strip("\n")
        elif line != "\n":
            coords=line.split(",")
            output='''%s\t%s\t%s\t%s\t%s\t%s''' % (currentChrom,currentGroup,coords[0],coords[1],coords[2],coords[3])
            o.write(output)
    c.close()
    o.close()

coordinatesToTSV(coordinateFile,outFile)
