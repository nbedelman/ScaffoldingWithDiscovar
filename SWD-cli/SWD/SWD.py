#!/usr/bin/env python

import argparse
import os
#
from comms import reOrderScaffolds
from comms import mafToBed as m2b
from comms import lastalOptions
from comms import lastdbOptions


parser = argparse.ArgumentParser(description='Scaffolding with DISCOVAR')

subparsers=parser.add_subparsers(help='sub-command help')

#genomeBuild subParser
parser_genomeBuild=subparsers.add_parser('genomeBuild',help='genomeBuild help')
parser_genomeBuild.add_argument('refGenome', type=str,
                    help='path to reference genome file')
parser_genomeBuild.add_argument('dbName', type=str,
                    help='name of database to create')
parser_genomeBuild.add_argument('config',nargs='?', type=str,
                    help='path to config file',default="lastdbConfig.py")
parser_genomeBuild.add_argument('--asChroms',action='store_true',
                    help='specifies that refGenome is a fasta file with one entry per chromosome. Default is that refGenome is a fasta file with one entry per scaffold.')
parser_genomeBuild.add_argument('--map', type=str,
                    help='path to map bed file')
parser_genomeBuild.set_defaults(which='genomeBuild')

#lastAlign subParser
parser_lastAlign=subparsers.add_parser('lastAlign',help='lastAlign help')
parser_lastAlign.add_argument('refLoc', type=str,
                    help='path to reference db, along with prefix. \
                    i.e. /path/to/reference/<dbName>')
parser_lastAlign.add_argument('query', type=str,
                    help='path to query sequence')
parser_lastAlign.add_argument('config',nargs='?', type=str,
                    help='path to config file',default="lastalConfig.py")
parser_lastAlign.set_defaults(which='lastAlign')

#reOrderScaffolds subParser
parser_reOrder=subparsers.add_parser('reOrder',help='reOrder help')
parser_reOrder.add_argument('bedDirectory', type=str,
                    help='path to directory containing bed files of contig \
                    alignments created with SWD mafToBed')
parser_reOrder.add_argument('map', type=str,
                    help='Original genome map file (AGP)')
parser_reOrder.add_argument('--refGenome',type=str,
                    help='Fasta file of reference genome')
parser_reOrder.add_argument('--newAssembly',type=str,
                    help='Fasta file of new assembly')
parser_reOrder.add_argument('--ungroupedChrom',default=None,type=str,
                    help='name of "chromosome" of ungrouped scaffolds, if it exists')
parser_reOrder.add_argument('--consecutiveOnly',action='store_true',
                    help='option to only join the obvious consecutive scaffolds')
parser_reOrder.add_argument('--combineMethod',default="first",type=str,
                    help='how to arrange the scaffolds when creating the new assembly')
parser_reOrder.add_argument('--reportDirectory',default=None,type=str,
                    help='option to only join the obvious consecutive scaffolds')
parser_reOrder.set_defaults(which='reOrder')


#mafToBed subParser
parser_mafToBed=subparsers.add_parser('mafToBed',help='mafToBed help')
parser_mafToBed.add_argument('alignment',type=str,
                            help='path to alignment in maf format')
parser_mafToBed.add_argument('assembly',type=str,
                             help='path to new assembly in fasta format')
parser_mafToBed.set_defaults(which='mafToBed')

#overviewToIntersects subParser
parser_overviewToIntersects=subparsers.add_parser('overviewToIntersects',help='overviewToIntersects help')
parser_overviewToIntersects.add_argument('overview',type=str,
                            help='path to overview bed file - should be \
                            something like <species>_allOverviews.bed')
parser_overviewToIntersects.add_argument('map',type=str,
                             help='path to map file in bed format')
parser_overviewToIntersects.add_argument('bedDirectory',type=str,
                             help='path to directory with all bed formats')
parser_overviewToIntersects.set_defaults(which='overviewToIntersects')

#Whichever parser was used, parse the arguments
arguments=parser.parse_args()

########### The genomeBuild Command ###########
def genomeBuild(refGenome,  dbName,config,asChroms,map=None):
    if not asChroms:
        chromGenome=os.path.basename(refGenome).split(".")[0]+"_chroms.fa"
        convertToChroms(map,refGenome, chromGenome)
        refGenome=chromGenome
    config=os.path.basename(config).split(".")[0]
    #print ("running lastdb with the following parameters:")
    #print ("config file:", config)
    options=lastdbOptions.getOptions(config)
    #print ("database name:", dbName)
    cmd='''sbatch genomeBuild.slurm "%s" %s %s''' % (options,dbName,refGenome)
    #print (cmd)
    os.system(cmd)
def convertToChroms(map,refGenome,chromGenome):
        print ("combining scaffolds into chromosomes")
        convertBed=os.path.basename(refGenome).split(".")[0]+"_tmpmap.bed"
        convertCMD='''awk -v OFS="\t" '{print $4,1,($3-$2)+1,$1,$5,$6}' %s > %s''' % (map,convertBed)
        print (convertCMD)
        os.system(convertCMD)
        combCMD='''bedtools getfasta -s -name -fi %s -bed %s -fo %s''' % (refGenome,convertBed,chromGenome+"_tmp")
        print (combCMD)
        os.system(combCMD)
        foldCMD='''fold -w 60 %s > %s ''' % (chromGenome+"_tmp", chromGenome+"_tmp2")
        os.system(foldCMD)
        #get rid of extra labels introduced in bedtools getfasta.
        chroms=open(chromGenome+"_tmp2","r")
        final=open(chromGenome,"w")
        used=[]
        for line in chroms:
            if '>' in line:
                if line not in used:
                    used.append(line)
                    final.write(line)
            else:
                final.write(line)
        chroms.close()
        final.close()
        #clean up the tmp files
        os.system('rm -f *tmp*')

if arguments.which=='genomeBuild':
    if not arguments.asChroms:
        genomeBuild(arguments.refGenome,  arguments.dbName,arguments.config,arguments.asChroms,arguments.map)
    else:
        genomeBuild(arguments.refGenome, arguments.dbName, arguments.config,arguments.asChroms)

###########the lastAlign Command#######

def lastAlign(config, refLoc, query):
    config=config.split("/")[-1].split(".")[0]
    print (config)
    options=lastalOptions.getOptions(config)
    cmd='''lastal %s %s %s''' % (options,refLoc,query)
    #cmd='''parallel-fastq "lastal %s %s " < %s ''' % (options,refLoc,query)
    print (cmd)
    os.system(cmd)

if arguments.which=='lastAlign':
    lastAlign(arguments.config, arguments.refLoc, arguments.query)

##############The mafToBed Command

def mafToBed(mafAlignment,discoOutput):
    m2b.runAll(mafAlignment,discoOutput)

if arguments.which=='mafToBed':
    mafToBed(arguments.alignment, arguments.assembly)

#############The overviewToIntersects Command#############
def O2I(bedFile,mafFile,bedDir):
    cmd='''overviewToIntersects.sh %s %s %s ''' % (bedFile,mafFile,bedDir)
    os.system(cmd)

if arguments.which=='overviewToIntersects':
    O2I(arguments.overview, arguments.map, arguments.bedDirectory)

#############The reOrderScaffolds Command#############

def reOrder(bedDirectory,map,refGenome,newAssembly,ungroupedChrom, combineMethod, reportDirectory, consecutiveOnly):
	reOrderScaffolds.runAll(bedDirectory,map,refGenome,newAssembly,ungroupedChrom, combineMethod, reportDirectory, consecutiveOnly)

if arguments.which=='reOrder':
	reOrder(arguments.bedDirectory,arguments.map,arguments.refGenome,arguments.newAssembly,arguments.ungroupedChrom, arguments.combineMethod, arguments.reportDirectory, arguments.consecutiveOnly)
