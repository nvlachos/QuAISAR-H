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
start = int(sys.argv[2])
end = int(sys.argv[3]) + 1
print(":"+str(start)+":"+str(end)+":")
print(main_record.id)
print(main_record.description)
search_DNA_seq = main_record.description[start:end]
print(search_DNA_seq)
