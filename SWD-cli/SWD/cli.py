#!/usr/bin/env python

import argparse
import os
from comms import hello
from comms import reOrderScaffolds
from comms import mafToBed
from comms import lastalOptions


parser = argparse.ArgumentParser(description='Scaffolding with DISCOVAR')

subparsers=parser.add_subparsers(help='sub-command help')

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
parser_reOrder.set_defaults(which='reOrder')


#mafToBed subParser
parser_mafToBed=subparsers.add_parser('mafToBed',help='mafToBed help')
parser_mafToBed.add_argument('alignment',type=str,
                            help='path to alignment in maf format')
parser_mafToBed.add_argument('assembly',type=str,
                             help='path to new assembly in fasta format')
parser_mafToBed.set_defaults(which='mafToBed')

#overviewToIntersects subParser
parser_mafToBed=subparsers.add_parser('overviewToIntersects',help='overviewToIntersects help')
parser_mafToBed.add_argument('overview',type=str,
                            help='path to overview bed file - should be \
                            something like <species>_allOverviews.bed')
parser_mafToBed.add_argument('map',type=str,
                             help='path to new assembly in fasta format')
parser_mafToBed.set_defaults(which='mafToBed')

#Whichever parser was used, parse the arguments
arguments=parser.parse_args()


###########the lastAlign Command#######

def lastAlign(config, refLoc, query):
    config=config.split("/")[-1].split(".")[0]
    print config
    options=lastalOptions.getOptions(config)
    cmd='''lastal %s %s %s''' % (options,refLoc,query)
    print cmd
    #os.system(cmd)

if arguments.which=='lastAlign':
    lastAlign(arguments.config, arguments.refLoc, arguments.query)

#############The reOrderScaffolds Command#############

def reOrder(bedDirectory,map,refGenome,newAssembly):
	reOrderScaffolds.runAll(bedDirectory,map,refGenome,newAssembly)

if arguments.which=='reOrder':
	reOrder(arguments.bedDirectory,arguments.map,arguments.refGenome,arguments.newAssembly)


##############The mafToBed Command

def mafToBed(mafAlignment,discoOutput):
    mafToBed.runAll(mafAlignment,discoOutput)

if arguments.which=='mafToBed':
    mafToBed(arguments.alignment, arguments.assembly)
