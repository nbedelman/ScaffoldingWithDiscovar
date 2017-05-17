import sys
from Bio import SeqIO

file=sys.argv[1]

fastaFile=open(file, "r")

used=[]

for line in fastaFile:
    if '>' in line:
	if line not in used:
	   used.append(line)
	   print line
    else:
	print line

fastaFile.close()


