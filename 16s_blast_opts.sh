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
. ./config.sh

# Load modules necessary for barrnap (that arent automatically loaded)
module load barrnap/0.8
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

#  Function to print out help blurb
show_help () {
	echo "Usage is ./16s_blast.sh -n sample_name -p run_id"
	echo "Output is saved to ${processed}/run_id/sample_name/16s"
}

options_found=0
while getopts ":h?n:p:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			project=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi


# Creates new output folder based on the universally set processed location from config.sh
OUTDATADIR=${processed}/${project}/${sample_name}/
if [ ! -d "${OUTDATADIR}/16s" ]; then
	echo "Creating $OUTDATADIR/16s"
	mkdir "${OUTDATADIR}/16s"
fi


# Function to create and add a "fasta" entry to the list of 16s hits
make_fasta() {
	header=">$3"
	rna_seq=""
	cstart=$4
	cstart=$(( cstart - 1 ))
	cstop=$5
	clength=$(( cstop - cstart + 1))
	match=0
	# Finds the matching contig and extracts sequence ##### REPLACE WITH SUBSEQUENCE ONCE IT CAN HANDLE MULTI_FASTAS !!!
	while IFS='' read -r line;
	do
		if [ "$line" == "${header}" ]; then
			match=1
			continue
		elif [[ "$line" = ">"* ]]; then
			match=0
		fi
		if [ $match -eq 1 ]; then
			rna_seq="$rna_seq$line"
		fi
	done < "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"
	# Extracts appropriate sequence from contig using start and stop positions
	rna="${rna_seq:$cstart:$clength}"
	# Adds new fasta entry to the file
	echo -e "${header}\n${rna_seq:$cstart:$clength}" >> ${processed}/${project}/${sample_name}/16s/${sample_name}_16s_rna_seqs.txt
}

owd=$(pwd)
cd ${OUTDATADIR}/16s

# Run barrnap to discover ribosomal sequences
barrnap --kingdom bac --threads ${procs} "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" > ${OUTDATADIR}/16s/${sample_name}_scaffolds_trimmed.fasta_rRNA_seqs.fasta

# Checks for successful output from barrnap, *rRNA_seqs.fasta
if [[ ! -s ${OUTDATADIR}/16s/${sample_name}_scaffolds_trimmed.fasta_rRNA_seqs.fasta ]]; then
	echo "rNA_seqs.fasta does NOT exist"
	exit 1
fi

# Checks barrnap output and finds all 16s hits and creates a fasta sequence to add to list of possible matches
lines=0
found_16s="false"
while IFS='' read -r line;
do
	if [ ${lines} -gt 0 ]; then
		contig=$(echo ${line} | cut -d' ' -f1)
		cstart=$(echo ${line} | cut -d' ' -f4)
		cstop=$(echo ${line} | cut -d' ' -f5)
		ribosome=$(echo ${line} | cut -d' ' -f9 | cut -d'=' -f3)
		if [ "${ribosome}" = "16S" ]; then
			# Replace with subsequence once it can handle multi-fastas
			#make_fasta $1 $2 $contig $cstart $cstop
			python3 ${shareScript}/get_subsequence.py -i "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" -s ${cstart} -e ${cstop} -t ${contig} >> ${processed}/${project}/${sample_name}/16s/${sample_name}_16s_rna_seqs.txt
			found_16s="true"
		fi
	fi
	lines=$((lines + 1))
done < "${OUTDATADIR}/16s/${sample_name}_scaffolds_trimmed.fasta_rRNA_seqs.fasta"

# Adds No hits found to output file in the case where no 16s ribosomal sequences were found
if [[ "${found_16s}" == "false" ]]; then
	echo -e "best_hit	${sample_name}	No_16s_sequences_found" > "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt"
	echo -e "largest_hit	${sample_name}	No_16s_sequences_found" >> "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt"
	exit
fi

# Blasts the NCBI database to find the closest hit to every entry in the 16s fasta list
###### MAX_TARGET_SEQS POSSIBLE ERROR
blastn -word_size 10 -task blastn -remote -db nt -max_hsps 1 -max_target_seqs 1 -query ${processed}/${project}/${sample_name}/16s/${sample_name}_16s_rna_seqs.txt -out ${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ssciname";
# Sorts the list based on sequence match length to find the largest hit
sort -k4 -n "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN" --reverse > "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN.sorted"

# Gets taxon info from the best (literal top) hit from the blast list
if [[ -s "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN" ]]; then
	me=$(whoami)
	accessions=$(head -n 1 "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN")
#	echo ${accessions}
	gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
	echo ${gb_acc}
	attempts=0
	# Will try getting info from entrez up to 5 times, as it has a higher chance of not finishing correctly on the first try
	while [[ ${attempts} -lt 5 ]]; do
		blast_id=$(python ${shareScript}/entrez_get_taxon_from_accession.py "${gb_acc}" "${me}@cdc.gov")
		if [[ ! -z ${blast_id} ]]; then
			break
		else
			attempts=$(( attempts + 1 ))
		fi
	done
	echo ${blast_id}
	if [[ -z ${blast_id} ]]; then
		blast_id="No_16s_matches_found"
	fi
	#blast_id=$(echo ${blast_id} | tr -d '\n')
	echo -e "best_hit	${sample_name}	${blast_id}" > "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt"
else
	echo "No remote blast file"
fi

# Gets taxon info from the largest hit from the blast list
if [[ -s "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN.sorted" ]]; then
	me=$(whoami)
	accessions=$(head -n 1 "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN.sorted")
	gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
	attempts=0
	# Will try getting info from entrez up to 5 times, as it has a higher chance of not finishing correctly on the first try
	while [[ ${attempts} -lt 5 ]]; do
		blast_id=$(python ${shareScript}/entrez_get_taxon_from_accession.py "${gb_acc}" "${me}@cdc.gov")
		if [[ ! -z ${blast_id} ]]; then
			break
		else
			attempts=$(( attempts + 1 ))
		fi
	done
	echo ${blast_id}
	if [[ -z ${blast_id} ]]; then
		blast_id="No_16s_matches_found"
	fi
#	blast_id$(echo ${blast_id} | tr -d '\n')
	echo -e "largest	${sample_name}	${blast_id}" >> "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt"
else
	echo "No sorted remote blast file"
fi

# Go back to original working directory
cd ${owd}

# Return modules to original versions
module unload perl/5.12.3
module load perl/5.22.1
module unload barrnap/0.8

#Script exited gracefully (unless something else inside failed)
exit 0
