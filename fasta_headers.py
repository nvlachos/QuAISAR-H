#!/usr/bin/env python3

'''
Changes fasta header from:
NODE_contig#_length_length#_cov_cov#
to
Name_contig#_length_length#


Usage: ./fasta_headers.py file.fasta output.fasta
'''

from Bio import SeqIO
import sys
import os

#print("Starting")

sequences = []
for record in SeqIO.parse(sys.argv[1],"fasta"):
    #print(record.id)
    name=os.path.basename(sys.argv[1]).split("_")[::-1]
    name=name[3:]
    name='_'.join(name[::-1])
    #print(name)
    record.id = record.id.split("_cov")[0].replace("NODE",name)
    #print(record.id)
    record.description = ""
#    print(record.description)
#    print(record)
    sequences.append(record)

SeqIO.write(sequences, sys.argv[2], "fasta")
