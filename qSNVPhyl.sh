#!/bin/sh -l

#$ -o qSNVPhyl.out
#$ -e qSNVPhyl.err
#$ -N qsnv
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# The wrapper script that runs the lowest possible run_SNVPhyl.sh script on the cluster
#
# Usage ./qSNVPhyl path_to_list_file output_directory project_identifier
#

# Checks for proper argumentation
if [[ $# -lt 3 ]]; then
	echo "Improper argument quantity supplied to $0, exiting"
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
	echo "No scheduling system, can not run qSNVPhyl.sh"
	#exit 1
fi
sleep 5

# Loop through and act on each sample name in the passed/provided list
counter=0
#echo "Starting counting"
while [[ -f ${shareScript}/run_SNVPhyl_${counter}.sh ]]; do
	counter=$(( counter + 1 ))
	echo ${counter}
	# Do we need to set a global variable for max SNVPhyls??? (instead of 10)
	if [[ ${counter} -ge 10 ]]; then
		#echo "Showing all jobs by MMB team..."
		qstat -u nvx4 -u nvd4 -u njr5 -u xku6 -u kqj9 -u nyx3 -u kbi5 -u yer1 -u vif0 -u hex1
		echo "Too many (10) SNVPhyls running already; please try again later (also check that there arent unused scripts in ${shareScript}/run_SNVPhyl_*.sh)"
		exit
	fi
done
#echo "Counted to ${counter}"

cp ${shareScript}/run_SNVPhyl_template.sh ${shareScript}/run_SNVPhyl_temp.sh
sed -i -e "s/run_SNVPhyl/run_SNVPhyl_${counter}/g" "${shareScript}/run_SNVPhyl_temp.sh"
sed -i -e "s/SNVPhyl_X/SNVPhyl_${counter}/g" "${shareScript}/run_SNVPhyl_temp.sh"
mv ${shareScript}/run_SNVPhyl_temp.sh ${shareScript}/run_SNVPhyl_${counter}.sh
echo "${shareScript}/run_SNVPhyl_${counter}.sh $@"
echo "Created and ran run_SNVPhyl_${counter}.sh"

qsub -sync y "${shareScript}/run_SNVPhyl_${counter}.sh" "$@"
rm ${shareScript}/run_SNVPhyl_${counter}.sh

global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
me=$(whoami)
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "qSNVPhyl.sh has completed running run_SNVPhyl_${counter}.sh" "${global_end_time}" | mail -s "SNVPhyl ${counter} analysis complete" "${submitter}@cdc.gov"

printf "%s %s" "qSNVPhyl.sh has completed running run_SNVPhyl_${counter}.sh " "${global_end_time}"
exit 0
