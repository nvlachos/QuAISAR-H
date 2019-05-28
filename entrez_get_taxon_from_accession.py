#!/usr/bin/env python
#  @ Author: Nick Vlachos
#  + Version .1
#  + 1 August 2017
#  Dependencies:  none

# Usage python ./entrez_get_taxon_from_accession.py accession_number your_email(for entrez tools)
#
# Requires module Entrez/E-utilities
#

from Bio import Entrez
import sys
import argparse

#Create an arg parser...someday
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Tool to retrieve taxonomy information from entrez using an accession number')
	parser.add_argument('-e', '--email', required=True, help='email of submitter, required by entrez')
	parser.add_argument('-a', '--accession', required=True, help='accession number to look up')
	return parser.parse_args()

args = parseArgs()
#Set the required email value to the supplied 2nd argument
Entrez.email = args.email
#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
handle = Entrez.efetch(db="nucleotide", id=sys.arg.accession, rettype="gb", retmode="text")
#Parses the returned output into lines
result=handle.read().split('\n')
#Goes through each line until it finds (and prints) the organism name that the accession number represents
for line in result:
	if 'ORGANISM' in line:
		print(' '.join(line.split()[1:]))
