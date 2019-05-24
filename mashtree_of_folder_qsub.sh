#!/bin/sh -l

#$ -o mashdist.out
#$ -e mashdist.err
#$ -N mashdist
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh

# Load necessary modules to run mashtree
module unload perl/5.22.1
module load perl/5.16.1-MT
module load mashtree/0.20

#
# Script to create mashtree of specified isolates that were processed by Quaisar pipeline
#
# Usage ./mashtree_of_list.sh path_to_assemblies output_filename extension_of_files_to_process
#

# create output directory if it does not exist
if [[ ! -d ${1} ]]; then
	mkdir -p ${1}
fi

# Call mashtree on all copied fasta
cd ${1}
mashtree.pl --numcpus ${procs} *.${3} --tempdir ${1}/temp > "${1}/${2}.dnd";

module unload perl/5.16.1-MT
module load perl/5.22.1

exit
