import sys
import Bio
import glob
import fileinput
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import unambiguous_dna
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from decimal import *
getcontext().prec = 3

def Gene_Counter(input_fasta):
    Gene_List = list(SeqIO.parse(input_fasta, 'fasta'))
    Gene_Count = len(Gene_List)
    return Gene_Count
    
def Duplicate_Gene_Remover(input_fasta, output_file):
    Gene_List = list(SeqIO.parse(input_fasta, 'fasta'))
    Gene_List_out = [Gene_List[0]]
    Gene_Adder = 0
    output_fasta = open(output_file, 'w')
    for genes in range(1, len(Gene_List)):
        for genes_2 in range(len(Gene_List_out)):
            if str(Gene_List[genes].seq) == str(Gene_List_out[genes_2].seq):
                Gene_Adder = 1
                break
            elif str(Gene_List[genes].seq.reverse_complement()) == str(Gene_List_out[genes_2].seq):
                Gene_Adder = 1
                break
            else:
                continue
        if Gene_Adder == 0:
##            print(Gene_List[genes].format('fasta'))
            Gene_List_out.append(Gene_List[genes])
        else:
            Gene_Adder = 0
    for output_genes in Gene_List_out:
        SeqIO.write(output_genes, output_fasta, 'fasta')
    output_fasta.close()

def Multi_Fasta_Combiner(index, new_file):
    """Reads in multiple fasta files and writes them out to a single new file"""
    File_List = glob.glob(index)
    output_handle = open(new_file, 'w') 
    for files in File_List:
        Gene_List = list(SeqIO.parse(files, 'fasta'))
        for genes in Gene_List:
            SeqIO.write(genes, output_handle, 'fasta')
    output_handle.close()

Multi_Fasta_Combiner(sys.argv[1] + '*.f*', sys.argv[2])
Duplicate_Gene_Remover(sys.argv[2], sys.argv[2])
