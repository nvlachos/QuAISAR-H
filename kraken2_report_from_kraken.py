import sys
import glob
import fileinput
import getpass
from Bio import Entrez

class taxon_Node:
	#initializing the variables
	name = ""
	count = 0
	taxID = 0
	parent = None
	taxLevel = "-"
	children = []

	#defining constructor
	def __init__(self, taxonName, taxonCount, taxonParent, taxonChildren, taxonID, taxonLevel):
		self.name = taxonName
		self.count = taxonCount
		self.parent = taxonParent
		self.children = taxonChildren
		self.taxonID = taxonID
		self.taxLevel = taxonLevel

	#defining constructor
	def __init__(self, taxonName, taxonCount, taxonParent, taxonID, taxonLevel):
		self.name = taxonName
		self.count = taxonCount
		self.parent = taxonParent
		self.taxonID = taxonID
		self.taxLevel = taxonLevel

	#defining class methods
	def showName(self):
		print(self.name)

	def showCount(self):
		print(self.count)

	def showParent(self):
		print(self.parent)

	def showtaxID(self):
		print(self.taxID)

	def showChildren(self):
		for child in self.children:
			print(child.showName)

	def addCounts(self, newReads):
		self.count += newReads

	def addChild(self, newChild):
		self.children.append(newChild)

	def Child_name(self, checkName):
		for child in self.children:
			if child.showName == checkName:
				return True

	def isChild_taxon(self, checkName):
		for child in self.children:
			if child.showtaxID == checkName:
				return True

	def find(self, checkName):
		for child in self.children:
			if child.showtaxID == checkName:
				return child
			else:
				child.find(checkName)
		return None

	#end of the class definition




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

# Script that will trim fasta files of any sequences that are smaller than the threshold
def get_mpa_string_From_NCBI(taxID):
	Entrez.email = getpass.getuser()
	#Creates the data structure from a pull from entrez nucleotide database using accession id with return type of genbank text mode
	handle = Entrez.efetch(db="taxonomy", id=taxID, mode="text", rettype="xml")
	#Parses the returned output into lines
	result= Entrez.read(handle)
	#Will have to change this region to be able to handle more descriptive lists downstream
	recognized_ranks={"superkingdom":"d", "kingdom":"k", "phylum":"p", "class":"c", "order":"o", "family":"f", "genus":"g", "species":"s"}
	for entry in result:
		#taxid = entry["Rank"]
		#print(entry)
		mpa_string=""
		for r in entry["LineageEx"]:
			#print(r)
			#print(r["Rank"])
			if r["Rank"] in recognized_ranks.keys():
				current_rank=recognized_ranks[r["Rank"]]
			#else:
				#current_rank="-"
				current_taxa=(r["ScientificName"])
				rank_and_taxa=current_rank+"__"+current_taxa
				#print(rank_and_taxa)
				mpa_string+=rank_and_taxa+"|"
		if entry["Rank"] in recognized_ranks.keys() or entry["Rank"] == "no rank":
			if entry["Rank"] == "no rank":
				current_rank="-"
			else:
				current_rank=recognized_ranks[entry["Rank"]]
		else:
			current_rank="-"
			current_taxa=(entry["ScientificName"])
			rank_and_taxa=current_rank+"__"+current_taxa
			#print(rank_and_taxa)
			mpa_string+=rank_and_taxa
		#else:
		#	print(entry["Rank"])
		return(mpa_string)

#def organize_mpas(input_kraken, output_mpa):
def organize_mpas():
	# kraken=open(input_kraken,'r')
	# line=kraken.readline().strip()
	# mpa_dict={}
	# mpa_counts={}
	# counter=0
	# while line != '':
	# 	line_sections = line.split("	")
	# 	contig_id = line_sections[1]
	# 	contig_taxID = line_sections[2]
	# 	if contig_taxID not in mpa_dict.keys():
	# 		print("Adding", contig_taxID)
	# 		mpa_dict[contig_taxID]=get_mpa_string_From_NCBI(contig_taxID)
	# 		mpa_counts[contig_taxID]=1
	# 	else:
	# 		if contig_taxID not in mpa_counts.keys():
	# 			#print("Doesnt exist???", contig_taxID)
	# 			mpa_counts[contig_taxID]=1
	# 		else:
	# 			#print("Incrementing:", contig_taxID)
	# 			mpa_counts[contig_taxID]+=1
	# 	line=kraken.readline().strip()
	# kraken.close()
	headNode = taxon_Node("unclassified", 0, None, 0, "u")
	headNode.addChild(taxon_Node("Test1", 23, taxon_Node("Test2", 34, None, 470, "G"), 469, "S"))
	print(headNode.find(469).showName)

	# mpa_taxon_counts={}
	# #print("mpa_dict length:", len(mpa_dict))
	# for key in mpa_dict.keys():
	# 	#print(key, mpa_dict[key])
	# 	taxons=mpa_dict[key].split("|")
	# 	if taxons is not None:
	# 		#print(taxons)
	# 		for i in range(0, len(taxons)):
	# 			#print(taxons[0:i+1])
	# 			if taxons[i] != "" and taxons[i][0:1] != "-":
	# 				#print(taxons[i])
	# 				if "|".join(taxons[0:i+1]) in mpa_taxon_counts:
	# 					#print("Incrementing", "|".join(taxons[0:i+1]), "from", mpa_taxon_counts["|".join(taxons[0:i+1])], "to",  mpa_taxon_counts["|".join(taxons[0:i+1])]+mpa_counts[key])
	# 					mpa_taxon_counts["|".join(taxons[0:i+1])]+=mpa_counts[key]
	# 				else:
	# 					#print("Creating", "|".join(taxons[0:i+1]), "at 1")
	# 					mpa_taxon_counts["|".join(taxons[0:i+1])]=mpa_counts[key]
	# 	else:
	# 		print("Taxons is none")
	#
	# print("mpa_taxon_counts length:", len(mpa_taxon_counts))
	# for key in sorted(mpa_taxon_counts.keys()):
	# 	print(key, mpa_taxon_counts[key])






	#print("mpa_counts length:", len(mpa_counts))
	#for key in mpa_counts.keys():
	#	print(key, mpa_counts[key])


#get_mpa_string_From_NCBI(470, blank_dick)

#organize_mpas(sys.argv[1], sys.argv[2])
organize_mpas()
