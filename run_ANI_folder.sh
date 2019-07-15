#!/bin/sh -l

#$ -o run_ANI_folder.out
#$ -e run_ANI_folder.err
#$ -N run_ANI_folder
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh
#. ${mod_changers}/pipeline_mods

ml Python3/3.5.2 pyani/0.2.7 mashtree/0.29

#
# Script to calculate the All vs. All average nucleotide identity of a folder of fasta files
#
# Usage ./run_ANI_folder.sh path_to_folder
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_ANI.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_ANI_folder.sh path_to_folder"
	exit 0
fi

python -V
echo "Running ALL vs ALL aniM on ${1} and placing results in ${1}/aniM"
#python "${shareScript}/pyani/average_nucleotide_identity.py" -i "${OUTDATADIR}/ANI/localANIDB" -o "${OUTDATADIR}/ANI/aniM" --write_excel
average_nucleotide_identity.py -i "${1}" -o "${1}/aniM"
