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
start = int(sys.argv[2])-1
end = int(sys.argv[3])
print(":"+str(start)+":"+str(end)+":")
print(main_record.id)
print(main_record.description)

search_DNA_seq = main_record.seq[start:end]
reverse_record = Seq(search_DNA_seq, generic_dna)
print(reverse_record.seq)
reverse_record=reverse_record.reverse_complement()
print(reverse_record.seq)

if len(sys.argv) > 4:
    if sys.argv[4] == "REVERSE" or sys.argv[4] == "R":
        search_DNA_seq=search_DNA_seq[::-1]
print(search_DNA_seq)
