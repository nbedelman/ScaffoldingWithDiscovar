#!/usr/bin/env python

#updated May 11, 2017
#updates include: score category 6 is not >50% in chrom instead of 75%
#Nate Edelman

from Bio import SeqIO
import sys
from operator import itemgetter

# mafAlignment=sys.argv[1]
# discoOutput=sys.argv[2]


def addToBedList(bedList, scaf, chromStart,chromEnd,id,length,strand,contig,conStart):
	'''takes a list of future bed entries, along with a number of attributes for a new entry.
	Adds the new bed attributes as a bed-ordered list'''
	bedList.append([scaf,chromStart-conStart,0,id,length,strand,chromStart,chromEnd,0])
	if length > bedList[0][4]:
		bedList[0]=[scaf,chromStart-conStart,0,contig,length,strand,chromStart,chromEnd,0]
	return bedList

class mafAttributes(object):
	def __init__(self,mafLine):
		line=mafLine.split()
		try:
			self.type=line[0]
			self.label=line[1]
			try:
				self.start=int(line[2])
			except ValueError:
				self.start=0
			try:
				self.end=self.start+int(line[3])
			except ValueError:
				self.end=0
			try:
				self.strand=line[4]
			except IndexError:
				self.strand="NA"
		except IndexError:
			self.type="NA"
			self.label="NA"
			self.start="NA"
			self.end="NA"
			self.length="NA"
			self.strand="NA"

	def getType(self):
		return self.type
	def getLabel(self):
		return self.label
	def getStart(self):
		return self.start
	def getEnd(self):
		return self.end
	def getLength(self):
		return (self.getEnd()-self.getStart())
	def getStrand(self):
		return self.strand


def mafToBedDict(mafFile):
	'''  Takes a maf file (like the output of LAST
	outputs a dictionary. Each key is a different query contig;  the value for each is
	a list of lists. Each single list contains ordered attributes to write to a bed file'''
	conDict={}
	current=[]
	id=0
	referenceStart=''
	with open(mafFile, "r") as file:
		#Skip all commented out lines -  there can be many at the beginning
		for line in file:
			if not "#" in line:
				maf=mafAttributes(line)
				if maf.getType()=='s':
					#Use the very first s line to define the reference pattern
					if referenceStart == '':
						referenceStart = maf.getLabel()[:3]
						#For every other line, we can now distinguish between reference and non-reference sequences
						#Use the reference line in a block to define the reference chromosome name and position
					if  referenceStart in maf.getLabel():
						#chrom=maf.getLabel()
						scaf=maf.getLabel()
						chromStart=maf.getStart()
						chromEnd=maf.getEnd()
						length=maf.getLength()
					#Use the query line to define the query name, position, and strand
					else:
						id += 1
						contig = maf.getLabel()
						if contig not in conDict.keys():
							conDict[contig] = [[0,0,0,0,0,0,0,0,0],]
							conStart= maf.getStart()
							#conEnd= maf.getEnd()
							strand = maf.getStrand()
							conDict[contig]=addToBedList(conDict[contig],scaf, chromStart,chromEnd,id,length,strand,contig,conStart)
						else:
							conStart= maf.getStart()
							strand = maf.getStrand()
							conDict[contig]=addToBedList(conDict[contig],scaf, chromStart,chromEnd,id,length,strand,contig,conStart)
	return conDict


def capBestAlign(fastaEntry,bedEntry):
    '''takes a specially-formatted bed-type entry and the fasta file that was used as its query
    returns a bed file where the entire query is placed relative to its best match Useful for visualizing alignments with IGV as well as properly re-ordering the scaffolds.'''
    bedEntry[2] = len(fastaEntry) + bedEntry[1]
    return bedEntry

def scoreAlign(record,entry):
    '''takes an alignment and the original query sequence.
    assigns a score based on the relative length of the alignment to the query'''
    length=entry[4]
    score=(float(length)/len(record))*1000
    entry[4]=score
    if score > 500:
        entry[8]="0,200,0"
    elif score > 250:
        entry[8]="128,255,0"
    elif score > 125:
        entry[8]="255,255,51"
    elif score > 67:
        entry[8]="255,153,51"
    else:
        entry[8]="255,0,0"
    return entry

def getBests(bedList):
    '''find and score the best alignment of a particular query sequence'''
    best=bedList[0]
    chrom=best[0]
    conStart=best[1]
    conEnd=best[2]
    conLength=abs(best[2]-best[1])
    insideCutoff=float(conLength)/20
    include=[best,]
    allInside=True
    sumChrom=0
    sumInside=0
    sumLength=0
    sort=sorted(bedList[1:], key=itemgetter(4), reverse=True)
    for entry in range(len(sort)):
        if sumLength<950:
            include.append(sort[entry])
            sumLength += sort[entry][4]
            if sort[entry][0] == chrom:
                sumChrom += sort[entry][4]
                if (abs((sort[entry][1] - conStart)) > insideCutoff) or ((sort[entry][2] - conEnd) > insideCutoff):
                    allInside = False
                else:
                    sumInside += sort[entry][4]
            else:
                allInside = False
        elif len(include) == 2:
            for repeat in range(entry,len(sort)):
                if sort[repeat][4] >= 900:
                    include.append(sort[repeat])
                    sumLength += sort[repeat][4]
                    allInside=False
                else:
                    break
            break
    category = categorizeBest(include,allInside, sumLength, sumInside, sumChrom)
    score=category[0]
    colorCode=category[1]
    include[0][4] = score
    include[0][8] = colorCode
    return include

def categorizeBest(include,allInside, sumLength, sumInside, sumChrom):
    if (len(include) == 2) and (sumInside >= 750):
        return [1000,"0,204,0"]
    elif allInside and (sumInside >= 750):
        return [900, "128,255,0"]
    elif sumLength >= 1900:
        if len(include) == 3:
            return [400, "51,255,255"]
        elif len(include) <= 5:
            return [300,"0,128,255"]
        else:
            return [200,"0,0,204"]
    elif sumInside >= 750:
        return [800,"255,255,0"]
    elif sumInside>=500:
        return [700,"255,153,51"]
    elif sumChrom >= 500:
        return [600, "255,102,102"]
    else:
        return [500, "204,0,0"]


def scoreAll(discoOutput,conDict):
    d = SeqIO.parse(discoOutput, "fasta")
    for record in d:
        if record.id in conDict.keys():
            conDict[record.id][0]=capBestAlign(record,conDict[record.id][0])
            for entry in range(1,len(conDict[record.id])):
                conDict[record.id][entry]=capBestAlign(record,conDict[record.id][entry])
                conDict[record.id][entry]=scoreAlign(record,conDict[record.id][entry])
            conDict[record.id]=getBests(conDict[record.id])
    return conDict

def writeToBed(bedList):
    name=bedList[0][3]
    f=open(name+".bed", "a+")
    for i in range(len(bedList)):
        toWrite=""
        for j in range(len(bedList[i])-1):
            toWrite+=(str(bedList[i][j])+" \t")
        toWrite+=(str(bedList[i][-1])+" \n")
        f.write(toWrite)
    f.close()


def runAll(mafAlignment,discoOutput):
    conDict = mafToBedDict(mafAlignment)
    scoreDict = scoreAll(discoOutput,conDict)
    for key in scoreDict.keys():
        writeToBed(scoreDict[key])

# runAll(mafAlignment,discoOutput)
