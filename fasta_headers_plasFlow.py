#!/usr/bin/env python3

'''
Changes fasta header from:
contig# length=length# depth=depthx
to
Name_contig#_length_length#_depth_depthx


Usage: ./fasta_headersplasFlow.py input.fasta output.fasta
'''

from Bio import SeqIO
import sys
import os

#print("Starting")

sequences = []
for record in SeqIO.parse(sys.argv[1],"fasta"):
    #print(record.id)
    name=os.path.basename(sys.argv[1]).split("_")[::-1]
    name=name[2:]
    name='_'.join(name[::-1])
    #print(name)
    #record.id = record.id.split("_cov")[0].replace("NODE",name)
    print(record)
    print(name)
    print(record.description.split(" ")[0])
    contig = record.description.split(" ")[0]
    record.description.split(" ")[1]
    length = record.description.split(" ")[1].split("=")[1]
    depth = record.description.split(" ")[2].split("=")[1]
    record.id = name+"_"+contig+"_length_"+length+"_depth_"+depth

    #print(record.id)
    record.description = ""
#    print(record.description)
#    print(record)
    sequences.append(record)

SeqIO.write(sequences, sys.argv[2], "fasta")
