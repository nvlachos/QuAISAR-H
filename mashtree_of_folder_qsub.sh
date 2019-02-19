#!/bin/sh -l

#$ -o mashdist.out
#$ -e mashdist.err
#$ -N mashdist
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh

# Load necessary modules to run mashtree
module unload perl/5.22.1
module load perl/5.16.1-MT

#
# Script to create mashtree of specified isolates that were processed by Quaisar pipeline
#
# Usage ./mashtree_of_list.sh path_to_assemblies output_filename extension_of_files_to_process
#

# create output directory if it does not exist
if [[ ! -d ${1} ]]; then
	mkdir -p ${1}
fi

# Copy over all fasta files from original locations to the output directory
while IFS= read -r line || [[ "$line" ]];  do
	sample_name=$(echo "${line}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${line}" | cut -d'/' -f1 | tr -d '[:space:]')
	cp ${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta ${2}
done < ${1}

# Call mashtree on all copied fasta
cd ${1}
mashtree.pl --numcpus ${procs} *.${3} --tempdir ${1}/temp > "${1}/${2}.dnd";

module unload perl/5.16.1-MT
module load perl/5.22.1

exit
