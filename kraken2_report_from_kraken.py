import sys
import glob
import fileinput
import getpass
from Bio import Entrez

class taxon_Node:
	#initializing the variables
	name = ""
	count = 0
	ID = 0
	parent = None
	level = "-"
	children = []

	#defining constructor
	def __init__(self, taxonName, taxonCount, taxonParent, taxonChildren, taxonID, taxonLevel):
		self.name = taxonName
		self.count = taxonCount
		self.parent = taxonParent
		if taxonChildren is not None:
			self.children = taxonChildren
		self.ID = taxonID
		self.level = taxonLevel

	#defining class methods
	def showName(self):
		#print("A")
		print(self.name)
		#print("B")

	def getName(self):
		return self.name

	def showCount(self):
		print(self.count)

	def getCount(self):
		return self.count

	def showParent(self):
		if self.parent is not None:
			print(self.parent)
		else:
			print("No parent assigned")

	def getParent(self):
		return self.parent

	def showtaxID(self):
		print(self.ID)

	def getTaxID(self):
		return taxID

	def showChildren(self):
		if self.children is None or len(self.children) == 0:
			print("No Children")
		else:
			print("Children:", len(self.children))
			for child in self.children:
				#print("Child:")
				child.showName()
				#print("dlihc")

	def getChildren(self):
		return self.children

	def addCounts(self, newReads):
		self.count += newReads

	def addChild(self, newChild):
		#print("Adding")
		#newChild.showName()
		#newChild.showtaxID()
		if self.getChildCount() == 0:
			#newChild.setParent(self)
			self.children=[newChild]
			#print("Children:", len(self.children), self.children[0].getName())
		else:
			#newChild.setParent(self)
			self.children.append(newChild)
			#print("Children:", len(self.children), self.children[0].getName(), self.children[len(self.children)-1].getName())
		#print("End Adding")

	def isChild_name(self, checkName):
		if self.children is not None:
			for child in self.children:
				if child.getName() == checkName:
					return True
		return False

	def isChild_taxon(self, checkName):
		if self.children is not None:
			for child in self.children:
				if child.getTaxID() == checkName:
					return True
		return False

	def find(self, checkName):
		if self.children is not None:
			for child in self.children:
				print(len(self.children), child.getTaxID(), checkName)
				if child.getTaxID() == checkName:
					return child
				else:
					self.find(checkName)
		return None

	def setParent(self, newParent):
		self.parent = newParent

	def getChildCount(self):
		if self.children is None or len(self.children) == 0:
			return 0
		else:
			return len(self.children)

	def showLevel(self):
		print(self.level)


	def print(self):
		print("Name:", self.name +"\nCounts:", str(self.count) +"\nID:", str(self.ID))
		print(self.getChildCount())
		self.showParent()
		self.showLevel()
		#print("End Print")
	#end of the class definition

#def organize_mpas(input_kraken, output_mpa):
def make_node_tree():
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
	headNode = taxon_Node("unclassified", 0, None, None, 0, "u")
	dNode = taxon_Node("Bacteria", 0, None, None, 2, "d")
	pNode = taxon_Node("Proteobacteria", 0, None, None, 1224, "p")
	cNode = taxon_Node("Gammaproteobacteria", 0, None, None, 1236, "c")
#	oNode = taxon_Node("Pseudomonadales", 0, None, 72274, "o")
	print("1-Show 3 Nodes Childrens")
	dNode.showChildren()
	dNode.showParent()
	dNode.print()
	print('2-Show headNode children')
	headNode.showChildren()
	headNode.showParent()
	headNode.print()
	print('3-Print dNode and link headnode to dNode and cNode')
	dNode.print()
	dNode.addChild(cNode)
	headNode.addChild(dNode)
	dNode.print()
	#headNode.addChild(dNode)
	print('4-Show headNode Children')
	headNode.showChildren()
	print('5-Show Parents of all children of headNode(', len(headNode.getChildren()), ")")
	for child in headNode.getChildren():
		if child.getParent() is not None:
			child.getParent().showName()
	print('6-Show grand children')
	for child in headNode.getChildren():
		child.showName()
		#print(child.getChildCount())
		child.showChildren()
	print('7')
	#headNode.showChildren()
	print('7')
	#headNode.find(2).showName()

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
make_node_tree()
