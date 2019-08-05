#!/bin/sh -l

#$ -o abl-template1.out
#$ -e abl-template1.err
#$ -N abl-template1
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
#. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_template.sh path_to_list path_for_output_file
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_template.sh path_for_list_file"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	while IFS= read -r line; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line:0:1}
		# Assign taxonomic level value from 4th value in line (1st-classification level, 2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "G" ]
		then
			genus=$(echo "${line}" | awk -F ' ' '{print $2}')
			# Only until ANI gets fixed
			if [[ ${genus} == "Clostridioides" ]]; then
				genus="Clostridium"
			fi
		fi
	done < "${processed}/${project}/${sample_name}/${sample_name}.tax"

	echo "${genus},${genus,,},${genus^}"

	if [[ -f "${processed}/${project}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus,,}.txt" ]]; then
		echo "Little exists"
		if [[ ! -f "${processed}/${project}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus^}.txt" ]]; then
			echo "Big does not"
			mv "${processed}/${project}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus,,}.txt" "${processed}/${project}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus^}.txt"
		fi
		rm "${processed}/${project}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus,,}.txt"
	fi
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list_template.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
