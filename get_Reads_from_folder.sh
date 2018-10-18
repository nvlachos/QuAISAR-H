#!/bin/sh -l

#$ -o get_Reads_from_Instruments.out
#$ -e get_Reads_from_Instruments.err
#$ -N get_Reads_from_Instruments
#$ -cwd
#$ -q all.q

# Import the config file with shortcuts and settings
. ./config.sh

#
# Will find all fastq.gz files within the given folder
#
# Usage ./get_Reads_from_folder.sh run_id folder_with_fastqs postfix_for_reads(1: _SX_RX_00X.fastq.gz 2: _RX.fastq.gz 3: _X.fastq.gz)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to get_Reads_from_folder.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty project name supplied to get_Reads_from_folder.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./get_Reads_from_folder.sh  run_id location_of_fastqs"
	echo "Output by default is downloaded to ${processed}/run_id and extracted to ${processed}/run_id/sample_name/FASTQs"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty source supplied to get_Reads_from_Folder.sh, exiting"
	exit 1
fi


# Sets folder to where files will be downloaded to
OUTDATADIR="${processed}/${1}"
if [ ! -d "${OUTDATADIR}" ]; then
	echo "Creating $OUTDATADIR"
	mkdir -p "${OUTDATADIR}"
fi
if [ -f "${OUTDATADIR}/${1}_list.txt" ]; then
	rm "${OUTDATADIR}/${1}_list.txt"
fi

####### Set trailing match pattern (everything after R in filename, to call unzip once for each pair) ######
match="${3}"

# Goes through given folder
echo "${2}"
for file in ${2}/*
do
	# Check if file is a zipped reads file
	if [[ "${file}" = *.gz ]]; then
		echo "filename: ${file}"
		# Gets full file name from path
		full_sample_name=${file##*/}
		# gets path from file
		source_path=$(dirname "${file}")
		# Extracts filename keeping only isolate ID, if it matches standard miseq naming
		if [[ "${match}" -eq 1 ]]; then
			short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f5- | rev)
			postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2,3 | rev)
			# Create an array out of the full sample name, delimited by _
			#IFS='_' read -r -a name_array <<< "${full_sample_name}"
			#long_name=${name_array[0]}_${name_array[1]}_${name_array[2]}_${name_array[3]}
		else
			short_name=$(echo "${full_sample_name}" | cut -d'_' -f1)
			postfix=$(echo "${full_sample_name}" | cut -d'_' -f2-)
		fi

		#long_name=$(echo "${full_sample_name}" | cut -d'_' -f1,2,3)
		echo "Short: ${short_name}"
		#echo "Does ${full_sample_name} match *${match}"

		# Skip file if it happens to be undetermined
		if [[ "${short_name}" == "Undetermined" ]]; then
			echo "found undetermined (${file})"
			continue
		# If the file matches the postfix given in the arguments proceed with moving and unzipping to the output directory
		else
			# Creates output folder
			if [ ! -d "${OUTDATADIR}/${short_name}" ]; then
				echo "Creating $OUTDATADIR/${short_name}"
				mkdir -p "${OUTDATADIR}/${short_name}"
				echo "Creating $OUTDATADIR/${short_name}/FASTQs"
				mkdir -p "${OUTDATADIR}/${short_name}/FASTQs"
			fi
			# Announces name of file being unzipped and then unzips it to the FASTQs folder for the matching sample name. Files are shortened to just name_R1_001.fastq or name_R2_001.fastq
			echo "Copying ${source_path}/${full_sample_name}"
			if [[ "${match}" -eq 1 ]] || [[ "${match}" -eq 2 ]]; then
				if [[ "${postfix}" = *"R1"* ]]; then
					echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					echo -e "${1}/${short_name}" >> "${OUTDATADIR}/${1}_list.txt"
				elif [[ "${postfix}" = *"R2"* ]]; then
					cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
				fi
			elif [[ "${match}" -eq 3 ]]; then
				if [[ "${postfix}" = *"1"* ]]; then
					echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					echo -e "${1}/${short_name}" >> "${OUTDATADIR}/${1}_list.txt"
				elif [[ "${postfix}" = *"2"* ]]; then
					cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
				fi
			else
				echo "Unrecognized postfix type, but how did it get this far?"
			fi
			# gunzip -c "${source_path}/${long_name}_R1_001.fastq.gz" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq"
			#cp "${source_path}/${long_name}_R1_001.fastq.gz" "${OUTDATADIR}/FASTQs/${short_name}_R1_001.fastq.gz"
			#echo "Unzipping ${long_name}_R2_001.fastq.gz"
			# gunzip -c "${source_path}/${long_name}_R2_001.fastq.gz" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq"
			#cp "${source_path}/${long_name}_R2_001.fastq.gz" "${OUTDATADIR}/FASTQs/${short_name}_R2_001.fastq.gz"
			# Add sample to list for current 'project' as it will be used by primary_processing to complete all downstream analyses
			#echo -e "${1}/${short_name}" >> "${OUTDATADIR}/${1}_list.txt"
		fi
	else
		echo "${file} is not a zipped read file, not acting on it"
	fi
	# Invert list so that the important isolates (for us at least) get run first
	if [[ -f "${OUTDATADIR}/${1}_list.txt" ]]; then
		sort -k2,2 -t'/' -r "${OUTDATADIR}/${1}_list.txt" -o "${OUTDATADIR}/${1}_list.txt"
	fi
done

#Script exited gracefully (unless something else inside failed)
exit 0
