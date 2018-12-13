#!/bin/sh -l

#$ -o run_ANI.out
#$ -e run_ANI.err
#$ -N run_ANI
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh
. ${mod_changers}/pipeline_mods

#
# Script to calculate the average nucleotide identity of a folder of fasta files
#
# Usage ./run_ANI.sh sample_name   path_to_folder_of_FASTAs
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_ANI.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_ANI.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_ANI.sh sample_name ani_database(which is also genus) species run_id list_of_samples_to_include(optional)"
	echo "Output is saved to in ${processed}/sample_name/ANI"
	exit 0
fi

python -V
echo "Running ALL vs ALL aniM on ${1} and placing results in ${1}/aniM"
#python "${shareScript}/pyani/average_nucleotide_identity.py" -i "${OUTDATADIR}/ANI/localANIDB" -o "${OUTDATADIR}/ANI/aniM" --write_excel
average_nucleotide_identity.py -i "${1}" -o "${1}/aniM"
