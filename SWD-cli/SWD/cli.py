#!/usr/bin/env python

import argparse
from comms import hello
from comms import reOrderScaffolds


parser = argparse.ArgumentParser(description='Scaffolding with DISCOVAR')

subparsers=parser.add_subparsers(help='sub-command help')
#Hello subParser
parser_hello=subparsers.add_parser('h',help='hello help')
parser_hello.add_argument('--Name', type=str,nargs='?',
                    help='Name of someone you\'d like to say hello to')
parser_hello.set_defaults(which='h')

#reOrderScaffolds subParser
parser_reOrder=subparsers.add_parser('reOrder',help='reOrder help')
parser_reOrder.add_argument('-b', '--bedDirectory',type=str,
                    help='path to directory containing bed files of contig alignments')
parser_reOrder.add_argument('-m', '--map',type=str,
                    help='Original genome map file (AGP)')
parser_reOrder.add_argument('-r', '--refGenome',type=str,
                    help='Fasta file of reference genome')
parser_reOrder.add_argument('-a', '--newAssembly',type=str,
                    help='Fasta file of new assembly')
parser_reOrder.set_defaults(which='reOrder')
arguments=parser.parse_args()


###########the Hello Command#######


def h(arg=''):
	hello.run(arg)

if arguments.which=='h':
	if arguments.Name:
		h(arguments.Name)
	else:
		h()


#############The reOrderScaffolds Command#############


def reOrder(bedDirectory,map,refGenome,newAssembly):
	reOrderScaffolds.runAll(bedDirectory,map,refGenome,newAssembly)

if arguments.which=='reOrder':
	reOrder(arguments.bedDirectory,arguments.map,arguments.refGenome,arguments.newAssembly)