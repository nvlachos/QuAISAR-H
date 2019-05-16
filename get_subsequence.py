#!/usr/bin/env python3

'''
Changes fasta header from:
contig# length=length# depth=depthx
to
Name_contig#_length_length#_depth_depthx


Usage: ./get_subsequence.py input.fasta output.fasta
'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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
reverse_record = SeqRecord(Seq(str(search_DNA_seq)), main_record.id, '', '')
#print(reverse_record.seq)
reverse_record=reverse_record.reverse_complement()
#print(reverse_record.seq)

if len(sys.argv) > 4:
    if len(sys.argv) == 7:
        start_sub = int(sys.argv[5])
        if start_sub > 0:
            start_sub-=1
        end_sub = int(sys.argv[6])
    else:
        start_sub = 0
        end_sub = len(str(reverse_record.seq))
    if sys.argv[4] == "REVERSE" or sys.argv[4] == "R":
        print(reverse_record.seq)[start_sub:end_sub]
        exit()
print(search_DNA_seq)
