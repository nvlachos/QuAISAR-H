#!/usr/bin/env python3

'''
Changes fasta header between:
NODE_contig#_length_length#_cov_cov#
and
Name_contig#_length_length#


Usage: ./fasta_headers.py -i file.fasta -o output.fasta [-r]
'''

from Bio import SeqIO
import sys
import os
import argparse

#print("Starting")
#Create an arg parser...someday
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to rename contigs in assemblies')
	parser.add_argument('-i', '--input', required=True, help='input fasta filename')
	parser.add_argument('-o', '--output', required=True, help='output filename')
	parser.add_argument('--reverse', help='returns formatted header to original', action='store_true')
	return parser.parse_args()

args=parseArgs()
sequences = []

if not args.reverse:
	#print("FORWARD")
	name=os.path.basename(args.input).split("_")[::-1]
	#print(name)
	name=name[3:]
	#print(name)
	name='_'.join(name[::-1])
	#print(name)
	for record in SeqIO.parse(args.input,"fasta"):
	    #print(record.id)
	    #print(name)
	    record.id = record.id.split("_cov")[0].replace("NODE",name)
	    #print(record.id)
	    record.description = ""
	#    print(record.description)
	#    print(record)
	    sequences.append(record)

	SeqIO.write(sequences, args.output, "fasta")
else:
	#print("REVERSE")
	name=os.path.basename(args.input).split("_")[::-1]
	#print(name)
	name=name[2:]
	#print(name)
	name='_'.join(name[::-1])
	#print(name)
	for record in SeqIO.parse(args.input,"fasta"):
	    #print(record.id)
	    #print(name)
	    record.id = record.id.replace(name,"NODE")+"_cov_X"
	    #print(record.id)
	    record.description = ""
		#print(record.description)
		#print(record)
	    sequences.append(record)

	SeqIO.write(sequences, args.output, "fasta")
