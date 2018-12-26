#!/bin/sh -l

#$ -o proj_alt_mlst.out
#$ -e proj_alt_mlst.err
#$ -N pam
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./run_csstar_proj_parser.sh list_name g/u (gapped/ungapped analysis ran) identity (80/95/98/99/100) output_name plasmid_identity(optional)
#
# Pulls out MLST, AR genes, and plasmid for the listed samples and consolidates them into one sheet (per category)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_csstar_proj_parser.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./project_alt_mlst.sh
	path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)
	output_location (will create folder with next parameter also if not existent)
	prefix_name_for_summary_file
	DB_to_look_for
	"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${share}"
	exit 0
elif [[ ! -f ${1} ]]; then
	echo "list does not exit...exiting"
	exit 1
elif [[ ! -d ${2}/${3} ]]; then
	mkdir -p "${2}/${3}"
fi

report_OUTDATADIR=${2}/${3}

# # Remove any pre-existing files from previous runs
if [[ ! -f ${report_OUTDATADIR}/${3}-alt_mlst_summary.txt ]]; then
	:
else
	rm ${report_OUTDATADIR}/${3}-alt-mlst_summary.txt
fi

# Loop through and act on each sample name in the passed/provided list
 while IFS= read -r line; do
 sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
 project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
 OUTDATADIR="${processed}/${project}/${sample_name}"


	# Pulls MLST type for sample and adds it to the summary file
	header=$(head -n 1 ${OUTDATADIR}/MLST/${sample_name}_${4}.mlst)
	IFS='	' read -r -a head_array <<< "$header"
	echo "head-${#head_array[@]}-${head_array[@]}"
	mlst=$(tail -n 1 ${OUTDATADIR}/MLST/${sample_name}_${4}.mlst)
	IFS='	' read -r -a mlst_array <<< "$mlst"
	echo "mlst-${#mlst_array[@]}-${mlst_array[@]}"

	#for index in ${#mlst_array}; do
	#	mlst_array[${index}]=$(echo "${mlst_array[${index}]}" | tr -d '[:space:]')
	#done

	echo -e "${mlst_array[1]}\t${project}\t${sample_name}\t${mlst_array[2]}\t${head_array[3]}(${mlst_array[3]})\t${head_array[4]}(${mlst_array[4]})\t${head_array[5]}(${mlst_array[5]})\t${head_array[6]}(${mlst_array[6]})\t${head_array[7]}(${mlst_array[7]})\t${head_array[8]}(${mlst_array[8]})\t${head_array[9]}(${mlst_array[9]})" >> ${report_OUTDATADIR}/${3}-alt-mlst_summary.txt

done < ${1}
