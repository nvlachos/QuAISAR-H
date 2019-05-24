#!/bin/sh -l

#$ -o qquaisar.out
#$ -e qquaisar.err
#$ -N qpp
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# The wrapper script that runs the lowest possible primary_processing.sh script on the cluster
#
# Usage ./primary_processing.sh
# -p/l/s(project/list/single) identifier: project pulls every sample from a particular "project" (or run identifier). Samples will be downloaded from instruments also.
# 							   			  list will only process samples found on the list
#							   			  sample will process the single sample supplied
#                                         *list and sample must have had the fastq files already downloaded and extracted
#

# Checks for proper argumentation
if [[ $# -lt 2 ]]; then
	echo "Improper argument quantity supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo -e "\\n\\n\\n"
			echo "Usage: ./qquaisar.sh [p/l/s] [project_name/list_of_samples.txt/sample_name] [n](optional)"
			echo "Must only be run on the cluster (ASPEN)"
			echo "n - no_download fastqs from instrument. Only apply to samples that have already had their FASTQ files placed into the proper locations"
			echo "p/l/s - project/list/single. What set of samples will the analysis be used on."
			echo "Must also give identifier to type of files, whether it is the project name, the text file containing a list of samples, or simply"
			echo "a single sample_id".
			echo "For samples in list or single mode - each name must include source project \"e.g. 170818_M02103_0075_000000000-B8VHY/1723664\""
			echo -e "\\n\\n\\n"
			exit 0
fi

# Not being run on cluster=no run
if [[ ${host} != "cluster"* ]]; then
	echo "No scheduling system, can not run qquaisar.sh"
	#exit 1
fi

sleep 5

# Loop through and act on each sample name in the passed/provided list
counter=0
#echo "Starting counting"
while [[ -f ${shareScript}/quaisar_${counter}.sh ]]; do
	counter=$(( counter + 1 ))
	echo ${counter}
	if [[ ${counter} -ge ${max_quaisars} ]]; then
		#echo "Showing all jobs by MMB team..."
		qstat -u nvx4 -u nvd4 -u njr5 -u xku6 -u kqj9 -u nyx3 -u kbi5 -u yer1 -u vif0 -u hex1
		echo "Too many quaisars running already; please try again later"
		exit
	fi
done
#echo "Counted to ${counter}"

cp ${shareScript}/quaisar_template.sh ${shareScript}/quaisar_${counter}.sh
sed -i -e "s/quaisar_/quaisar_${counter}/g" "${shareScript}/quaisar_${counter}.sh"
sed -i -e "s/quasX/quas${counter}/g" "${shareScript}/quaisar_${counter}.sh"
#mv ${shareScript}/primary_processing_temp.sh ${shareScript}/primary_processing_${counter}.sh
echo "${shareScript}/quaisar_${counter}.sh $@"
echo "Created and ran quaisar_${counter}.sh"

qsub -sync y "${shareScript}/quaisar_${counter}.sh" "$@"
rm ${shareScript}/quaisar_${counter}.sh

global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "qquasingle.sh has completed running quaisar_${counter}.sh " "${global_end_time}"
exit 0
