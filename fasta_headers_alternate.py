#!/usr/bin/env python3

#
# Description: Changes headers in SPAdes assembly fasta from contig# length=length# depth=depthx to Name_contig#_length_length#_depth_depthx
#
# Usage: ./fasta_headers_plasFlow.py -i input.fasta -o output.fasta
#
# Output location: parameter
#
# Modules required: Biopython must be available in python instance
#
# v1.0 (10/3/2019)
#
# Created by Erisa Sula (nvd4@cdc.gov)
#

from Bio import SeqIO
import sys
import os
import argparse

#print("Starting")
#Create an arg parser...someday
def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to rename contigs in NCBI assemblies')
    parser.add_argument('-i', '--input', required=True, help='input fasta filename')
    parser.add_argument('-o', '--output', required=True, help='output filename')
    return parser.parse_args()

args=parseArgs()
sequences = []


#print("FORWARD")
name=os.path.basename(args.input).split('_')[0]

print(name)
#name=name[3:]
#print(name)
#name='-'.join(name[::])
print("Name=", name)
for record in SeqIO.parse(args.input,"fasta"):
    #print(record.id)
    #print(record.name)
    #print(name)
    #record.id = record.id.split("_cov")[0].replace("NODE",name)
    entrails = record.id.split("|")   #[6].split("_")[-1]
    #print(*entrails, sep = "\n")
    contig_num=entrails[2].split("_")[5]
    seq_length=len(record)
    #print("Contig=", contig_num)
    print(name+"_"+str(contig_num)+"_length_"+str(seq_length)+"|"+record.id)
    record.id = name+"_"+str(contig_num)+"_length_"+str(seq_length)+"|"+record.id
    #print(record.id)
    record.description = ""
#    print(record.description)
#    print(record)
    sequences.append(record)
SeqIO.write(sequences, args.output, "fasta")
