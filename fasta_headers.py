#!/usr/bin/env python3

'''
Changes fasta header from:
NODE_contig#_length_length#_cov_cov#
to
Name_contig#_length_length#


Usage: ./fasta_headers.py -i file.fasta -o output.fasta
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
	return parser.parse_args()

args=parseArgs()
sequences = []
for record in SeqIO.parse(args.input,"fasta"):
    #print(record.id)
    name=os.path.basename(args.input).split("_")[::-1]
    name=name[3:]
    name='_'.join(name[::-1])
    #print(name)
    record.id = record.id.split("_cov")[0].replace("NODE",name)
    #print(record.id)
    record.description = ""
#    print(record.description)
#    print(record)
    sequences.append(record)

SeqIO.write(sequences, args.output, "fasta")
