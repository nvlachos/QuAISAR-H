#!/usr/bin/env python3

'''
Changes fasta header from:
contig# length=length# depth=depthx
to
Name_contig#_length_length#_depth_depthx


Usage: ./get_subsequence.py input.fasta output.fasta
'''

from Bio import SeqIO
import sys
import os

#print("Starting")

sequence = ""
main_record = SeqIO.read(sys.argv[1],"fasta")
print(":"+str(sys.argv[2])+":"+str(sys.argv[3])+":")
search_DNA_seq = str(main_record)[sys.argv[2]:sys.argv[3]]
print(search_DNA_seq)
