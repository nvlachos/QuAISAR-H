import sys
import glob
import fileinput
import getpass


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


def translate(input_kraken, output_labels):
	kraken=open(input_kraken,'r')
	line=kraken.readline().strip()
	tax_tree_dict={}
	label_lines=[]
	while line != '':
		line_sections = line.split("	")
		contig_id = line_sections[1]
		contig_taxID = line_sections[2]
		if contig_taxID in tax_tree_dict.keys():
			label_lines.append(contig_id+"	"+tax_tree_dict[contig_taxID])
		else:
			tax_tree_dict[contig_taxID]=get_Taxon_Tree_From_NCBI(contig_taxID)
		label_lines.append(contig_id+"	"+tax_tree_dict[contig_taxID])
	print("Lines:", len(label_lines))

translate(sys.argv[1])
