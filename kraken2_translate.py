import sys
import glob
import fileinput
import getpass
from Bio import Entrez

# Script that will trim fasta files of any sequences that are smaller than the threshold
def get_Taxon_Tree_From_NCBI(taxID):
	Entrez.email = getpass.getuser()
	#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
	handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text", rettype="xml")
	#Parses the returned output into lines
	result= Entrez.read(handle)
	#Goes through each line until it finds (and prints) the organism name that the accession number represents
	for taxon in result:
		taxid = taxon["TaxId"]
		name = taxon["ScientificName"]
		lineage=["root"]
		for t in taxon["LineageEx"]:
			lineage.append(t["ScientificName"])
		lineage.append(name)
		#print("%s\t|\t%s\t|\t%s" % (taxid, name, ";".join(lineage)))
		return ";".join(lineage)
		#print(' '.join(line.split()[1:]))

# Script that will trim fasta files of any sequences that are smaller than the threshold
def get_mpa_string_From_NCBI(taxID, mpa_dict):
	Entrez.email = getpass.getuser()
	#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
	handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text", rettype="xml")
	#Parses the returned output into lines
	result= Entrez.read(handle)
	#Goes through each line until it finds (and prints) the organism name that the accession number represents
	recognized_ranks={"superkingdom":"d", "kingdom":"k", "phylum":"p", "class":"c", "order":"o", "family":"f", "genus":"g", "species":"s", "species_group":"x"}
	for entry in result:
		#taxid = entry["Rank"]
		mpa_string=""
		for r in entry["LineageEx"]:
			print(r)
			print(r["Rank"])
			if r["Rank"] in recognized_ranks.keys():
				current_rank=recognized_ranks[r["Rank"]]
			#else:
				#current_rank="-"
				current_taxa=(r["ScientificName"])
				rank_and_taxa=current_rank+"__"+current_taxa
				print(rank_and_taxa)
				mpa_string+=rank_and_taxa+"|"
		if entry["Rank"] in recognized_ranks.keys():
			current_rank=recognized_ranks[entry["Rank"]]
		#else:
			#current_rank="-"
			current_taxa=(entry["ScientificName"])
			rank_and_taxa=current_rank+"__"+current_taxa
			print(rank_and_taxa)
			mpa_string+=rank_and_taxa
		print(mpa_string)

def translate(input_kraken, output_labels):
	kraken=open(input_kraken,'r')
	line=kraken.readline().strip()
	tax_tree_dict={}
	mpa_dict={}
	label_lines=[]
	counter=0
	while line != '':
		line_sections = line.split("	")
		contig_id = line_sections[1]
		contig_taxID = line_sections[2]
		if contig_taxID not in tax_tree_dict.keys():
			tax_tree_dict[contig_taxID]=get_Taxon_Tree_From_NCBI(contig_taxID)
		print(str(counter)+":"+contig_id+"	"+tax_tree_dict[contig_taxID])
		label_lines.append(contig_id+"	"+tax_tree_dict[contig_taxID])
		line=kraken.readline().strip()
		counter+=1
	kraken.close()
	label_file=open(output_labels, 'w')
	label_file.write("\n".join(label_lines))
	print("Lines:", len(label_lines))
	for line in label_lines:
		print(line)

#translate(sys.argv[1], sys.argv[2])

blank_dick={}
get_mpa_string_From_NCBI(470, blank_dick)
