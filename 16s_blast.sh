#!/bin/sh -l

#$ -o 16s_blast.out
#$ -e 16s_blast.err
#$ -N 16s_blast
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi

#. ./module_changers/list_modules.sh

# Load modules necessary for barrnap (that arent automatically loaded)
module unload perl/5.22.1
module load perl/5.12.3

#
# Creates a species prediction based on blasting the largest and also best hit of the suggested 16s sequences found using barrnap
# Usage ./16s_blast.sh   sample_name   run_id
#
# Required modules: barrnap/0.8
# Sub-required modules (loaded by required modules): hmmer/3.1b2
#
# Required modules that arent automatically loaded: perl 5.12.3 (not 5.22.1)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to 16s_blast.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to 16s_blast.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./16s_blast.sh   sample_name   run_id"
	echo "Output is saved to ${processed}/run_id/sample_name/16s"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty run_id supplied to 16s_blast.sh, exiting"
	exit 1
fi

# Creates new output folder based on the universally set processed location from config.sh
OUTDATADIR=${processed}/${2}/${1}
if [ ! -d "${OUTDATADIR}/16s" ]; then
	echo "Creating $OUTDATADIR/16s"
	mkdir "${OUTDATADIR}/16s"
fi


# Function to create and add a "fasta" entry to the list of 16s hits
make_fasta() {
#	echo "Made it in fasta maker - $1,$2,$3,$4,$5"
	header=">$3"
	rna_seq=""
	cstart=$4
	cstart=$(( cstart - 1 ))
	cstop=$5
	clength=$(( cstop - cstart + 1))
	match=0
	# Finds the matching contig and extracts sequence
	while IFS='' read -r line;
	do
#		echo "test before-${line}=${header}?"
		if [ "$line" == "${header}" ]; then
			match=1
			continue
		elif [[ "$line" = ">"* ]]; then
			match=0
		fi
		if [ $match -eq 1 ]; then
		#echo "${line}"
		#echo "-${rna_seq}"
			rna_seq="$rna_seq$line"
		fi
	done < "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta"
#	echo "$cstart:$clength"
	# Extracts appropriate sequence from contig using start and stop positions
	rna="${rna_seq:$cstart:$clength}"
#	echo "short_seq=$rna"
	# Adds new fasta entry to the file
	echo -e "${header}\n${rna_seq:$cstart:$clength}" >> ${processed}/${2}/${1}/16s/${1}_16s_rna_seqs.txt
}

owd=$(pwd)
cd ${OUTDATADIR}/16s

#Run rnammer to predict 16s RNA sequence NO LONGER USED, but keeping for posterity
#rnammer -S bac -m ssu -xml ${OUTDATADIR}/16s/${1}_scaffolds_trimmed.fasta_rRNA_seqs.xml -gff ${OUTDATADIR}/16s/${1}_scaffolds_trimmed.fasta_rRNA_seqs.gff -h ${OUTDATADIR}/16s/${1}_scaffolds_trimmed.fasta_rRNA_seqs.hmmreport -f ${OUTDATADIR}/16s/${1}_scaffolds_trimmed.fasta_rRNA_seqs.fasta < ${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta

# Run barrnap to discover ribosomal sequences
barrnap --kingdom bac --threads ${procs} "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" > ${OUTDATADIR}/16s/${1}_scaffolds_trimmed.fasta_rRNA_seqs.fasta

# Checks for successful output from barrnap, *rRNA_seqs.fasta
if [[ ! -s ${OUTDATADIR}/16s/${1}_scaffolds_trimmed.fasta_rRNA_seqs.fasta ]]; then
	echo "rNA_seqs.fasta does NOT exist"
	exit 1
fi

# Checks brrnap output and finds all 16s hits and creates a fasta sequence to add to list of possible hits
lines=0
found_16s="false"
while IFS='' read -r line;
do
	if [ ${lines} -gt 0 ]; then
		contig=$(echo ${line} | cut -d' ' -f1)
		cstart=$(echo ${line} | cut -d' ' -f4)
		cstop=$(echo ${line} | cut -d' ' -f5)
		ribosome=$(echo ${line} | cut -d' ' -f9 | cut -d'=' -f3)
#		echo "ribo-$ribosome"
		if [ "${ribosome}" = "16S" ]; then
			make_fasta $1 $2 $contig $cstart $cstop
			found_16s="true"
		fi
	fi
	lines=$((lines + 1))
done < "${OUTDATADIR}/16s/${1}_scaffolds_trimmed.fasta_rRNA_seqs.fasta"

# Adds No hits found to output file in the case where no 16s ribosomal sequences were found
if [[ "${found_16s}" == "false" ]]; then
	echo -e "best_hit	${1}	No_16s_sequences_found" > "${OUTDATADIR}/16s/${1}_16s_blast_id.txt"
	echo -e "largest_hit	${1}	No_16s_sequences_found" >> "${OUTDATADIR}/16s/${1}_16s_blast_id.txt"
	exit
fi

# Blasts the NCBI database to find the closest hit to every entry in the 16s fasta list
###### MAX_TARGET_SEQS POSSIBLE ERROR
blastn -word_size 10 -task blastn -remote -db nt -max_hsps 1 -max_target_seqs 1 -query ${processed}/${2}/${1}/16s/${1}_16s_rna_seqs.txt -out ${OUTDATADIR}/16s/${1}.nt.RemoteBLASTN -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ssciname";
# Sorts the list based on sequence match length to find the largest hit
sort -k4 -n "${OUTDATADIR}/16s/${1}.nt.RemoteBLASTN" --reverse > "${OUTDATADIR}/16s/${1}.nt.RemoteBLASTN.sorted"

# Gets taxon info from the best (literal top) hit from the blast list
if [[ -s "${OUTDATADIR}/16s/${1}.nt.RemoteBLASTN" ]]; then
	me=$(whoami)
	accessions=$(head -n 1 "${OUTDATADIR}/16s/${1}.nt.RemoteBLASTN")
#	echo ${accessions}
	gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
	echo ${gb_acc}
	blast_id=$(python ${shareScript}/entrez_get_taxon_from_accession.py "${gb_acc}" "${me}@cdc.gov")
	echo ${blast_id}
	if [[ -z ${blast_id} ]]; then
		blast_id="No_16s_matches_found"
	fi
	#blast_id=$(echo ${blast_id} | tr -d '\n')
	echo -e "best_hit	${1}	${blast_id}" > "${OUTDATADIR}/16s/${1}_16s_blast_id.txt"
else
	echo "No remote blast file"
fi

# Gets taxon info from the largest hit from the blast list
if [[ -s "${OUTDATADIR}/16s/${1}.nt.RemoteBLASTN.sorted" ]]; then
	me=$(whoami)
	accessions=$(head -n 1 "${OUTDATADIR}/16s/${1}.nt.RemoteBLASTN.sorted")
#	echo ${accessions}
	gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
	blast_id=$(python ${shareScript}/entrez_get_taxon_from_accession.py "${gb_acc}" "${me}@cdc.gov")
	echo ${blast_id}
	if [[ -z ${blast_id} ]]; then
		blast_id="No_16s_matches_found"
	fi
#	blast_id$(echo ${blast_id} | tr -d '\n')
	echo -e "largest	${1}	${blast_id}" >> "${OUTDATADIR}/16s/${1}_16s_blast_id.txt"
else
	echo "No sorted remote blast file"
fi

# Go back to original working directory
cd ${owd}

# Return modules to original versions
module load perl/5.12.3
module load perl/5.22.1

#Script exited gracefully (unless something else inside failed)
exit 0
