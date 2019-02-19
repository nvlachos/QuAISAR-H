#!/bin/sh -l

#$ -o act_by_list_template.out
#$ -e act_by_list_template.err
#$ -N abl-template
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list_AR_completion_check.sh path_to_list ResGANNOT_identifier(YYYYMMDD)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_AR_completion_check.sh path_to_list_file ResGANNOT_Identifier(YYYYMMDD)"
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list
counter=1
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	"${shareScript}/determine_taxID.sh" "${sample_name}" "${project}"
	while IFS= read -r line;
	do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			species=$(echo "${line}" | awk -F '	' '{print $2}')
		elif [ "${first}" = "G" ]
		then
			genus=$(echo "${line}" | awk -F ' ' '{print $2}')
		fi
	done < "${processed}/${project}/${sample_name}/${sample_name}.tax"

	"${shareScript}/run_MLST.sh" ${sample_name} ${project}
	if [[ "${genus}_${species}" = "Acinetobacter_baumannii" ]]; then
		"${shareScript}/run_MLST.sh" "${sample_name}" "${project}" "-f" "abaumannii"
	elif [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
		"${shareScript}/run_MLST.sh" "${sample_name}" "${project}" "-f" "ecoli_2"
	fi

	"${shareScript}/validate_piperun.sh" ${sample_name} ${project} > "${processed}/${project}/${sample_name}/"
	echo "${counter}:${sample_name}/${project} completed"
	counter=$(( counter + 1 ))
done < "${1}"

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
printf "%s %s" "Act_by_list.sh has completed updating MLST and pipeline_stats " "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
