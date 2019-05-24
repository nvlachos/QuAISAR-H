#!/bin/sh -l

#$ -o run_SPAdes.out
#$ -e run_SPAdes.err
#$ -N run_SPAdes
#$ -cwd
#$ -q all.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Runs SPAdes on sample to align reads into best possible assembly
#
# Usage ./run_SCCmecFinder.sh sample_name	run_id
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_SPAdes.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_SCCmecFinder.sh sample_name	run_id"
	echo "Output by default is sent to ${processed}/miseq_run_id/sample_name/SCCmec"
	exit 0
elif [ -z "${2}" ]; then
	echo "Empty project id supplied to run_SPAdes.sh, exiting"
	exit 1
fi

if [[ ! -d "${processed}/${2}/${1}/SCCmec" ]]; then
	mkdir -p "${processed}/${2}/${1}/SCCmec"
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}/SCCmec"

python2 "${shareScript}/SCCmecFinder_v4.py"

#Script exited gracefully (unless something else inside failed)
exit 0
