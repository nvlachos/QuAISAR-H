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
	echo "Improper argument quantity supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo -e "\\n\\n\\n"
			echo "Usage: ./qSNVPhylsh path_to_list output_directory project_identifier"
			echo "Must only be run on the cluster (ASPEN)"
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
while [[ -f ${shareScript}/run_SNVPhyl_${counter}.sh ]]; do
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

cp ${shareScript}/run_SNVPhyl_template.sh ${shareScript}/run_SNVPhyl_temp.sh
sed -i -e "s/run_SNVPhyl/run_SNVPhyl_${counter}/g" "${shareScript}/run_SNVPhyl_temp.sh"
mv ${shareScript}/run_SNVPhyl_temp.sh ${shareScript}/run_SNVPhyl_${counter}.sh
echo "${shareScript}/run_SNVPhyl_${counter}.sh $@"
echo "Created and ran run_SNVPhyl_${counter}.sh"

qsub -sync y "${shareScript}/run_SNVPhyl_${counter}.sh" "$@"
rm ${shareScript}/run_SNVPhyl_${counter}.sh

global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "qSNVPhyl.sh has completed running run_SNVPhyl_${counter}.sh " "${global_end_time}"
exit 0
