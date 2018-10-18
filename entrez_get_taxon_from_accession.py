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

#Create an arg parser...someday

#Set the required email value to the supplied 2nd argument
Entrez.email = sys.argv[2]
#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
handle = Entrez.efetch(db="nucleotide", id=sys.argv[1], rettype="gb", retmode="text")
#Parses the returned output into lines
result=handle.read().split('\n')
#Goes through each line until it finds (and prints) the organism name that the accession number represents
for line in result:
	if 'ORGANISM' in line:
		print(' '.join(line.split()[1:]))
