#!/bin/sh -l

#$ -o gss.out
#$ -e gss.err
#$ -N gss
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: A small script to get overall status of sample and specifically assembly status
#
# Usage: ./get_sample_status.sh sample_name project
#
# Output location: standard out
#
# Modules required: None
#
# v1.0 (10/28/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Shows a brief uasge/help section if -h option used as first argument
if [[ "$1" = "-h" ]]; then
	echo -e "\\n\\n\\n"
	echo "Usage: ./get_sample_status.sh sample_name project"
	echo -e "\\n\\n\\n"
	exit 0
# Checks for proper argumentation
elif [[ $# -gt 2 ]]; then
	echo "Improper argument quantity supplied to $0, exiting"
	exit 1
fi

if [[ ! -f "${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt" ]]; then
	echo "Pipeline stats file does not exist, can not continue"
else
	assembly="Unknown"
	completion="Unknown"
	while IFS= read -r var; do
		if [[ "${var}" == "Assembly            :"* ]];
			assembly=$(echo "${var}" | cut -d':' -f2 | tr -d [:space:])
			break
		fi
	done < "${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt"
	completion=$(tail -n1 "${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt" | cut -d' ' -f5)
	echo "${completion}:${assembly}"
fi

#Script exited gracefully (unless something else inside failed)
exit 0
