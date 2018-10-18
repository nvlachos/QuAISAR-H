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
#. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/module_changers/pipeline_mods
#. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/module_changers/list_modules.sh

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
	echo "Usage is ./run_plasmidFinder.sh  sample_name run_id output_folder(either plasmid or plasmid_on_plasmidAssembly) (-i number_minimum_identity, optional) (-f to force against all databases, optional)"
	echo "Output by default is ${processed}/miseq_run_id/sample_name/plasmid"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty run_id supplied to run_plasFlow.sh, exiting"
	exit 1
fi

module load PlasFlow/1.1
module load Python/3.5.4
. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/module_changers/list_modules.sh

if [[ ! -d "${processed}/${2}/${1}/plasFlow" ]]; then
	mkdir "${processed}/${2}/${1}/plasFlow"
fi

if [[ -s "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta" ]]; then
	perl "${shareScript}/removeShortContigs_2000.pl" "${processed}/${project}/${1}/Assembly/${1}_scaffolds_trimmed.fasta"
	mv "${processed}/${project}/${1}/Assembly/${1}_scaffolds_trimmed.fasta.TRIMMED.fasta" "${processed}/${project}/${1}/plasFlow/${1}_scaffolds_trimmed_2000.fasta"
	PlasFlow.py --input "${processed}/${2}/${1}/plasFlow/${1}_scaffolds_trimmed_2000.fasta" --output "${processed}/${2}/${1}/plasFlow/${1}_plasFlow.tsv"
fi
