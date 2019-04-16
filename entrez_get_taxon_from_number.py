#!/usr/bin/env python
#  @ Author: Nick Vlachos
#  + Version .1
#  + 1 August 2017
#  Dependencies:  none

# Usage python ./entrez_get_taxon_from_number.py ncbi-taxonomy-number your_email(for entrez tools)
#
# Requires module Entrez/E-utilities
#

from Bio import Entrez
import sys

#Create an arg parser...someday

#Set the required email value to the supplied 2nd argument
Entrez.email = sys.argv[2]
#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
handle = Entrez.efetch(db="taxonomy", id=sys.argv[1], mode="text", rettype="xml")
#Parses the returned output into lines
result= Entrez.read(handle)
#Goes through each line until it finds (and prints) the organism name that the accession number represents
for taxon in result:
	taxid = taxon["TaxId"]
	name = taxon["ScientificName"]
	lineage=[]
	for t in taxon["LineageEx"]:
		lineage.append(t["ScientificName"])
	lineage.append(taxid)
	print("%s\t|\t%s\t|\t%s" % (taxid, name, ";".join(lineage)))
	#print(' '.join(line.split()[1:]))
