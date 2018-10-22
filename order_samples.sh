#!/bin/sh -l

#$ -o order_samples.out
#$ -e order_samples.err
#$ -N order_samples
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
# No MODS needed

#
# Creates a txt list file that contains run and samples for a given MiSeq run that matches the order of the output for the MMB_Seq log
#
# Usage ./order_samples.sh MiSeq_Run_ID
#
# Output will be /MiseqAnalysisFiles/Run_ID/Run_ID_list_ordered.txt
#


#### #copy MMBSEQLog to local
> "${processed}/${1}/${1}_list_ordered.txt"

echo "${processed}/${1}/2019_MMBSeq_Log.xlsx"

# Copy the newest log file to the local directory
cp "${shareScript}/2019_MMBSeq_Log.xlsx" "${processed}/${1}/2019_MMBSeq_Log.xlsx"

# Convert log file to csv format for searchability
`soffice --headless --convert-to csv "${processed}/${1}/2019_MMBSeq_Log.xlsx" --outdir "${processed}/${1}/" > "${processed}/${1}/convert.out"`
echo "Excel file: 2019_MMBSeq_Log.xlsx has been converted to CSV"

# Parse log file csv until run_if matches
while IFS= read -r var; do
	# Check the format of the city/state column in the log file to determine how many tabs need to be used to find run_id in line
	#echo "checking ${var}"
	city_state=$(echo "${var}" | cut -d',' -f16)
	#echo "|^| "${#city_state}" : "${city_state}
	if [[ ${#city_state} -eq 2 ]] || [[ ${city_state} = "unknown" ]] || [[ ${city_state} = "Unknown" ]] || [[ ${city_state} = "UNKNOWN" ]] || [[ "${city_state}" = "" ]]; then
		#echo "Doing if"
		line_project=$(echo "${var}" | cut -d',' -f20)
	else
		#echo "doing else"
		line_project=$(echo "${var}" | cut -d',' -f21)
	fi
	#echo "${line_project}:${1}"
	# If the run_id matches, then add ID to list (automatically placing them in the proper order)
	if [[ "${line_project}" = "${1}" ]]; then
		line_id=$(echo "${var}" | cut -d',' -f3)
		#echo "${1}/${line_id}"
		echo "${1}/${line_id}" >> "${processed}/${1}/${1}_list_ordered.txt"
	else
		#echo "Not in ${1}"
		:
	fi
done < ${processed}/${1}/2019_MMBSeq_Log.csv

# Remove intermediate files from sorting
rm -r ${processed}/${1}/sorted_summaries.txt
rm -r ${processed}/${1}/2019_MMBSeq_Log.csv
rm -r ${processed}/${1}/2019_MMBSeq_Log.xlsx
rm -r ${processed}/${1}/convert.out

# Check if the sorted file has content, else delete it since something went wrong
if [[ ! -s "${processed}/${1}/${1}_list_ordered.txt" ]]; then
	echo "Isolates were not able to be sorted, something wrong with MiSeq Log entries or list file, or....?"
	rm -r ${processed}/${1}/${1}_list_ordered.txt
	exit
else
	echo "sorted file contains entries"
	:
fi
