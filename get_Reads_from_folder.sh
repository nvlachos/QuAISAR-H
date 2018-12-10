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
# Will find all fastq.gz files within the given folder
#
# Usage ./get_Reads_from_folder.sh run_id folder_with_fastqs postfix_for_reads(1:_l001_SX_RX_00X.fastq.gz 2: _RX.fastq.gz 3: _X.fastq.gz 4: _RX_00X.fastq.gz)
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
	if [[ "${file}" = *.gz ]] || [[ "${file}" = *.fastq ]]; then
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
			if [[ "${match}" -eq 1 ]] || [[ "${match}" -eq 4 ]] || ; then
				if [[ "${postfix}" = *"R1_001.fast"* ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"

					fi
					echo -e "${1}/${short_name}" >> "${OUTDATADIR}/${1}_list.txt"
				elif [[ "${postfix}" = *"R2_001.fast"* ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
						cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
						gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
					fi
				fi
			elif [[ "${match}" -eq 3 ]]; then
				if [[ "${postfix}" = *"1.fast"* ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					fi
					echo -e "${1}/${short_name}" >> "${OUTDATADIR}/${1}_list.txt"
				elif [[ "${postfix}" = *"2.fast"* ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
						cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
						gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
					fi
				fi
			elif [[ "${match}" -eq 2 ]]; then
				if [[ "${postfix}" = *"R1.fast"* ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
					fi
					echo -e "${1}/${short_name}" >> "${OUTDATADIR}/${1}_list.txt"
				elif [[ "${postfix}" = *"R2.fast"* ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
						cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
						gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz"
						clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001_unclumped.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
					fi
				fi
			else
				echo "Unrecognized postfix type, but how did it get this far?"
			fi
		fi
	else
		echo "${file} is not a zipped read file, not acting on it"
	fi
done

# Invert list so that the important isolates (for us at least) get run first
if [[ -f "${OUTDATADIR}/${1}_list.txt" ]]; then
	sort -k2,2 -t'/' -r "${OUTDATADIR}/${1}_list.txt" -o "${OUTDATADIR}/${1}_list.txt"
fi

#Script exited gracefully (unless something else inside failed)
exit 0
