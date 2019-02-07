#!/bin/sh -l

#$ -o get_Reads_from_Instruments.out
#$ -e get_Reads_from_Instruments.err
#$ -N get_Reads_from_Instruments
#$ -cwd
#$ -q all.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Will find all fastq.gz files matching run ID on any of our machines (listed in config.sh)
#
# Usage ./get_Reads_from_Instruments.sh project_name
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to get_Reads_from_Instrument.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty project name supplied to get_Reads_from_Instrument.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./get_Reads_from_Instrument.sh  miseq_run_id"
	echo "Output by default is extracted to ${processed}/miseq_run_id/sample_name/FASTQs"
	exit 0
fi

project_found="false"

# Checks to see if a sample list file exists in the destination. Deletes it if found because a new one will be made from what is downloaded
OUTDATADIR="${processed}/${1}"
if [[ -s "${processed}/${1}/${1}_list.txt" ]]; then
	rm "${processed}/${1}/${1}_list.txt"
fi

# Goes through each folder of the machines listed in the config file (Currently 3 miseqs)
for instrument in "${all_instruments[@]}"
do
	# Goes through each subfolder of the current instrument
	for run_folder in "${instrument}"/*
	do
		# Gets folder names in current directory
		run_id=${run_folder##*/}
		#echo "${run_id} - ${1}"
		# If folder name matches project name
		if [[ "${run_id}" = "${1}" ]]; then
			# Print that a match was found
			echo "Found project ${run_id} in ${instrument}"
			# Go through every file in the Basecalls folder of the found folder (all files will only be fastq.gzs)
			for file in ${run_folder}/Data/Intensities/BaseCalls/*.fastq.gz
			do
				# Drops path info and keeps only filename
				full_sample_name=${file##*/}
				source_path=$(dirname "${file}")
				# Short name will contain only the isolate ID for a sample
				short_name=$(echo "${full_sample_name}" | cut -d'_' -f1)
				# Make an array out of the elements in a miseq run id, delimited by underscores
				IFS='_' read -r -a name_array <<< "${full_sample_name}"
				# Long name will contain the isolate ID as well as the S# from the run
				long_name=${name_array[0]}_${name_array[1]}
				# Check if sample id is undetermined so that these are not downloaded
				if [[ "${short_name}" == "Undetermined" ]]; then
					continue
				# Acts upon the R1 of each sample (as long as there is a matching R2) to create a sample and FASTQs folder and then download reads into
				elif [[ ${full_sample_name} == *"L001_R1_001.fastq.gz" ]] && [[ -f "${source_path}/${long_name}_L001_R2_001.fastq.gz" ]]; then
					# Creates destination folder
					if [ ! -d "${OUTDATADIR}/${short_name}" ]; then
						mkdir -p "${OUTDATADIR}/${short_name}"
						mkdir -p "${OUTDATADIR}/${short_name}/FASTQs"
					fi
					#clumpify.sh in="${source_path}/${long_name}_L001_R1_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					#clumpify.sh in="${source_path}/${long_name}_L001_R2_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
					clumpify.sh in1="${source_path}/${long_name}_L001_R1_001.fastq.gz" in2="${source_path}/${long_name}_L001_R2_001.fastq.gz" out1="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out2="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" reorder
					gunzip -c "${source_path}/${long_name}_L001_R1_001.fastq.gz" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq"
					gunzip -c "${source_path}/${long_name}_L001_R2_001.fastq.gz" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq"
					# Add sample to list for current 'project' as it will be used by the pipeline to complete all downstream analyses
					echo -e "${1}/${short_name}" >> "${processed}/${1}/${1}_list.txt"
				# If only the R2 is present, then create sample and FASTQs folder to download into...also notify maintenance that there is a problem
				elif [[ ${full_sample_name} == *"L001_R2_001.fastq.gz" ]] && [[ ! -f "${source_path}/${long_name}_L001_R1_001.fastq.gz" ]]; then
					# Creates destination folder
					if [ ! -d "${OUTDATADIR}/${short_name}" ]; then
						echo "Creating $OUTDATADIR/${short_name}"
						mkdir -p "${OUTDATADIR}/${short_name}"
						mkdir -p "${OUTDATADIR}/${short_name}/FASTQs"
					fi
					# Print a warning tag to alert user about missing read
					echo "




					R2 ONLY - Unzipping ${long_name}_L001_R2_001.fastq.gz




					"
					clumpify.sh in="${source_path}/${long_name}_L001_R2_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
					gunzip -c "${source_path}/${long_name}_L001_R2_001.fastq.gz" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq"
					# Add sample to list for current 'project' as it will be used by primary_processing to complete all downstream analyses
					echo -e "${1}/${short_name} - R2 ONLY from ${file}" >> "${processed}/${1}/${1}_list.txt"
					echo -e "${1}/${short_name} - R2 ONLY" >> "${shareScript}/maintenance_To_Do.txt"
				# If only the R1 is present, then create sample and FASTQs folder to download into...also notify maintenance that there is a problem
				elif [[ ${full_sample_name} == *"L001_R1_001.fastq.gz" ]] && [[ ! -f "${source_path}/${long_name}_L001_R2_001.fastq.gz" ]]; then
					# Creates destination folder
					if [ ! -d "${OUTDATADIR}/${short_name}" ]; then
						echo "Creating $OUTDATADIR/${short_name}"
						mkdir -p "${OUTDATADIR}/${short_name}"
						mkdir -p "${OUTDATADIR}/${short_name}/FASTQs"
					fi
					# Print a warning tag to alert user about missing read
					echo "



					R1 ONLY - Unzipping ${long_name}_L001_R1_001.fastq.gz



					"
					clumpify.sh in="${source_path}/${long_name}_L001_R1_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					gunzip -c "${source_path}/${long_name}_L001_R1_001.fastq.gz" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq"
					# Add sample to list for current 'project' as it will be used by primary_processing to complete all downstream analyses
					echo -e "${1}/${short_name} - R1 ONLY" >> "${processed}/${1}/${1}_list.txt"
					echo -e "${1}/${short_name} - R1 ONLY from ${file}" >> "${shareScript}/maintenance_To_Do.txt"
				elif [[ ${full_sample_name} == *"L001_I1_001.fastq.gz" ]] || [[ ${full_sample_name} == *"L001_I2_001.fastq.gz" ]]; then
					echo "Index reads are not used currently: (${full_sample_name}"
					continue
				fi
			done
			project_found="true"
		fi
	done
	# Invert list so that the important isolates (for us at least) get run first
	if [[ -f "${OUTDATADIR}/${1}_list.txt" ]]; then
		sort -k2,2 -t'/' -r "${OUTDATADIR}/${1}_list.txt" -o "${OUTDATADIR}/${1}_list.txt"
	fi
done

# Report if a matching run_id was found in the available instruments
if [[ "${project_found}" = "false" ]]; then
	echo "Project not found in ${all_instruments[@]}"
fi
#Script exited gracefully (unless something else inside failed)
exit 0
