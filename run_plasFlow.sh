#!/bin/sh -l

#$ -o get_run_plasFlow.out
#$ -e get_run_plasFlow.err
#$ -N get_run_plasFlow
#$ -cwd
#$ -q all.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#. "${mod_changers}/pipeline_mods"
#. ./module_changers/pipeline_mods
#. ./module_changers/list_modules.sh

#
# Will attempt to find any plasmids in sample
#
# Usage ./run_plasFlow.sh sample_name run_id
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_plasmidFinder.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_plasmidFinder.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_plasmFlow.sh  sample_name run_id"
	echo "Output by default is ${processed}/miseq_run_id/sample_name/plasmFlow"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty run_id supplied to run_plasFlow.sh, exiting"
	exit 1
fi

module load PlasFlow/1.1
module load Python/3.5.4
. ./module_changers/list_modules.sh

if [[ ! -d "${processed}/${2}/${1}/plasFlow" ]]; then
	mkdir "${processed}/${2}/${1}/plasFlow"
fi

if [[ -s "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta" ]]; then
	python2 "${shareScript}/removeShortContigs.py" "${processed}/${project}/${1}/Assembly/scaffolds.fasta" "2000"
	${shareScript}/fasta_headers.py "${processed}/${project}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${processed}/${project}/${1}/plasFlow/${1}_scaffolds_trimmed_2000.fasta"
	rm -r "${processed}/${project}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta"
	PlasFlow.py --input "${processed}/${2}/${1}/plasFlow/${1}_scaffolds_trimmed_2000.fasta" --output "${processed}/${2}/${1}/plasFlow/${1}_plasFlow.tsv"
fi
