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
import argparse
import sys
import os

#print("Starting")


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='input file', required=True, dest='input_file')
parser.add_argument('-s', '--start', help='start position', required=True, type=int, dest='start')
parser.add_argument('-e', '--end', help='end position', required=True, type=int, dest='end')
parser.add_argument('-t', '--title', help='header', required=False, type=str, dest='title')
parser.add_argument('-r', '--reverse', help='reverse complement strand?', required=False, type=bool, dest='reverse', default=False)
parser.add_argument('-b', '--begin', help='rev-comp begin position', required=False, type=int, dest='begin', default=0)
parser.add_argument('-f', '--finish', help='rev-comp finish position', required=False, type=int, dest='finish', default=0)
parameters=parser.parse_args()

with open(parameters.input_file, 'r') as myfile:
    data=myfile.read().replace('\n', ' ')
print(data.split(' ').count(">"))
exit()

match_record=None

for record in SeqIO.parse(parameters.input_file, "fasta"):
    if record.id.starts_with(parameters.header):
            match_record = record
            print("Found matching header")

if match_record == None:
    print("No matching header found, exiting...")

sequence = ""
#match_record = SeqIO.read(sys.argv[1],"fasta")
match_record = SeqIO.read(parameters.input_file, "fasta")
#start = int(sys.argv[2])-1
#end = int(sys.argv[3])
print(":"+str(parameters.start)+":"+str(parameters.end)+":")
print(match_record.id)
print(match_record.description)

search_DNA_seq = match_record.seq[parameters.start:parameters.end]
reverse_record = SeqRecord(Seq(str(search_DNA_seq)), match_record.id, '', '')
#print(reverse_record.seq)
reverse_record=reverse_record.reverse_complement()
#print(reverse_record.seq)

if parameters.reverse:
    start_sub = parameters.begin
    if start_sub > 0:
        start_sub-=1
    if parameters.finish <= parameters.begin:
        print("Reverse finish is less than reverse begin???")
        exit()
    else:
        end_sub = parameters.finish
    print(reverse_record.seq)[start_sub:end_sub]
    exit()
print(search_DNA_seq)
