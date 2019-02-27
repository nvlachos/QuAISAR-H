#!/bin/bash -l

#$ -o parqua.out
#$ -e parqua.err
#$ -N parqua
#$ -cwd
#$ -q all.q

# Copy config file into working directory to allow changes to made to output directory if necessary
shareScript=$(pwd)
echo "${shareScript}"
config_counter=0
while true
do
	if [[ "${config_counter}" -gt 9 ]]; then
		echo "Already 10 parallel quaisar sets running, please wait until one finishes...exiting"
		exit 324
	fi
	if [[ ! -f "${shareScript}/config_${config_counter}.sh" ]]; then
		if [[ -f "${shareScript}/config_template.sh" ]]; then
			echo "Trying to copy config_template.sh to config_${config_counter}"
			cp "${shareScript}/config_template.sh" "${shareScript}/config_${config_counter}.sh"
			break
		else
			echo "config_template.sh does not exist, cannot copy and must exit..."
			exit 333
		fi
	else
		config_counter=$(( config_counter + 1 ))
	fi
done

#Import the config file with shortcuts and settings
. ${shareScript}/config_${config_counter}.sh

#Import the module file that loads all necessary mods
. ${mod_changers}/pipeline_mods





#
# The wrapper script that runs all the tools that have been designated as necessary (and some others that are typically run also)
#
# Usage ./primary_processing.sh
# -i path (1,2,3)						: Path to directory containing all zipped fastq files, must include numerical identifier of postfix naming scheme (1=_SX_RX_00X_fastq.gz. 2=_RX.fastq.gz. 3=X.fastq.gz)
# -p/l/s(project/list/single) identifier: project pulls every sample from a particular "project" (or run identifier). Samples will be downloaded from instruments also.
# 							   			  list will only process samples found on the list
#							   			  sample will process the single sample supplied
#                                         *list and sample must have had the fastq files already downloaded and extracted
#

#Print out which type of machine the script is running on (Biolinux or Aspen as an interactive session or node based)
if [ "${host}" = "biolinux" ];
then
	echo "Running pipeline on Biolinux"
elif [ "${host}" = "aspen_login" ];
then
	echo "Running pipeline on Aspen interactive node"
elif [[ "${host}" = "cluster"* ]];
then
	echo "Running pipeline on Aspen ${host}"
fi

# Remove all modules as to not interfere when the new list is load in
yes | bash "${mod_changers}/clear_mods.sh"

#load all proper modules and output currently loaded ones
echo "Loading necessary modules for pipeline"
. ${mod_changers}/pipeline_mods
${mod_changers}/list_modules.sh

# Checking BASH version
if [ "${BASH_VERSINFO}" -lt 4 ];
then
	echo "Sorry, you need at least bash-4.0 to run this script." >&2;
	exit 1;
fi

# Checking for proper number of arguments from command line
if [[ $# -lt 1  || $# -gt 7 ]]; then
        echo "Usage: ./primary_processing.sh [i] [o] [d/n] [p/l/s] [project_name/list_of_samples.txt/sample_name]"
        echo "You have used $# args"
        exit 3
fi

# Checks the arguments (more to come)
nopts=$#
do_download=false
global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
requestor=$(whoami)
PROJECT="${requestor}_${global_time}"
BASEDIR="${processed}"

for ((i=1 ; i <= nopts ; i++)); do
	#echo "${1} ${2}"
	case "${1}" in
		#Help/Usage section
		-h | --help)
			echo -e "\\n\\n\\n"
			echo "Usage: ./parallel_quisar.sh [i] directory number [o] directory name [p] project_name"
			echo "i - input folder location containing all ZIPPED FASTQ files, must also include additional parameter after this indicating what postfix files will be"
			echo "1,2,3 - options that correlate to postfix files types (1=_SX_RX_00X_fastq.gz. 2=_RX.fastq.gz. 3=X.fastq.gz)"
			echo "o - output folder where to put project folder"
			echo "n/d - no_download/download fastqs from instrument. Only applies to project analysis (Samples can not be downloaded yet in list or single mode)"
			echo "p/l/s - project/list/single. What set of samples will the analysis be used on."
			echo "Must also give identifier to type of files, whether it is the project name, the text file containing a list of samples, or simply"
			echo "a single sample_id".
			echo "For samples in list or single mode - each name must include source project \"e.g. 170818_M02103_0075_000000000-B8VHY/1723664\""
			echo -e "\\n\\n\\n"
			exit 0
			;;
		#Gets name of folder that FASTA files will be in
		-i | --in-dir)
			INDATADIR="$2"
			if [[ -d  ${INDATADIR} ]]; then
				do_download="true"
				list_path="${BASEDIR}/${PROJECT}/${PROJECT}_list.txt"
			else
					echo "${INDATADIR} does not exist...exiting"
					exit 1
			fi
			indir_set="true"
			postfix="$3"
			#is_full_run="false"
			#echo "$INDATADIR $2"
			shift 3
			;;
		#Gets output directory name of folder that all output files will be stored
		-o | --out-dir)
			BASEDIR="$2"
			PROJECT="$3"
			shift 3

			echo "processed=${BASEDIR}" >> "${shareScript}/config.sh"
			. ${shareScript}/config.sh
			echo "${processed}"
			list_path="${BASEDIR}/${PROJECT}/${PROJECT}_list.txt"
			if [[ ! -d ${BASEDIR} ]]; then
				mkdir -p ${BASEDIR}
			fi
			;;
		#Checks for (project) name of folder that all output files will be stored
		-p | --project-name)
			PROJECT="$2"
			is_proj="true"
			#Copies over and unzips all 'Determined' zipped FASTQ files from the sequencing instrument for the requested project_id
			for instrument in "${all_instruments[@]}"
			do
				# Goes through each subfolder of the current instrument
				for run_folder in "${instrument}"/*
				do
					# Gets folder names in current directory
					run_id=${run_folder##*/}
					#echo "${run_id} - ${PROJECT}"
					# If folder name matches project name
					if [[ "${run_id}" = "${PROJECT}" ]]; then
						# Print that a match was found
						echo "Found project ${run_id} in ${instrument}"
						# Go through every file in the Basecalls folder of the found folder (all files will only be fastq.gzs)
						INDATADIR="${instrument}/${run_id}/Data/Intensities/BaseCalls"
						break
					fi
				done
			done
			if [[ ! -d "${INDATADIR}" ]]; then
				echo "No FOLDER ${INDATADIR} exists on any MiSeq instrument, exiting"
				exit 123
			fi
			if [[ ${BASEDIR} = "${requestor}_${global_time}" ]]; then
				BASEDIR="${processed}"
			fi
			list_path="${processed}/${PROJECT}/${PROJECT}_list.txt"
			postfix=1
			shift 2
			;;
		# #Tells the script to run analyses from already downloaded fastq files
		# -n | --no_download)
		# 	do_download=false
		# 	shift
		# 	;;
		# #Tells the script that only the files found in the attached list need to be run
		# -l | --list)
		# 	list_path="$2"
		# 	quick_list=$(echo "${2}" | cut -d'.' -f1)
		# 	#INDATADIR="${processed}/${PROJECT}" # NOT USED yet in list mode
		# 	if [[ -z ${BASEDIR} ]]; then
		# 		BASEDIR="${processed}"
		# 	fi
		# 	do_download="false"
		# 	PROJECT="list_${quick_list}"
		# 	is_full_run="false"
		# 	shift 2
		# 	;;
		# #Tells the script that only the single isolate needs to be run
		# -s | --single)
		# 	echo "${2}/${3}" > "./tempList_${global_time}.txt"
		# 	list_path="./tempList_${global_time}.txt"
		# 	INDATADIR="${processed}/${3}"
		# 	if [[ -z ${BASEDIR} ]]; then
		# 		BASEDIR="${processed}"
		# 	fi
			#trn=$(echo "${2}" | sed 's/\//-/')
			#PROJECT="single_${trn}"
			# PROJECT="run_for_sample_${3}"
			# do_download="false"
			# is_full_run="false"
			# shift 3
			# ;;
		#Captures any other characters in the args
		\?)
			echo "ERROR: ${BOLD}$2${NORM} is not a valid argument" >&2
			usage
			exit 1
			;;
	esac
done

# Short print out summary of run settings
echo -e "Source folder: ${INDATADIR}\\nOutput folder: ${BASEDIR}\\nList based analysis:  ${list_path}"

# Checks that a full FASTQ source path is given
if [ "${INDATADIR:0:1}" != "/" ] && [ "${INDATADIR:0:1}" != "$" ]; then
	echo "${INDATADIR}"
	echo 'ERROR: The full path was not specified.' >&2
	exit 1
fi

# Copies reads from source location to working directory and creates a list of IDs
"${shareScript}/get_Reads_from_folder.sh" "${PROJECT}" "${INDATADIR}" "${postfix}"

# Loops through list file to create an array of all isolates to run through pipeline
declare -a file_list=()
while IFS= read -r file;
do
	echo "Found: ${file}"
	file=$(echo "${file}" | tr -d '[:space:]')
	file_list+=("${file}")
done < "${list_path}"

# Displays number and names of files found to analyze
if [[ ${#file_list[@]} -gt 1 ]]; then
	echo "Will analyze these ${#file_list[@]} files: ${file_list[*]}"
elif [[ ${#file_list[@]} -eq 1 ]]; then
	echo "Will analyze this file: ${file_list[0]}"
else
	echo "No files found in ${list_path}"
fi

run_start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

#Each file in the list is checked individually for successful completion and added then added to the log for the run
mkdir -p "${Quaisar_H_log_directory}/${PROJECT}_on_${run_start_time}"
log_dir="${Quaisar_H_log_directory}/${PROJECT}_on_${run_start_time}"

#Get the time the run started to use as the identifier
outarray=()
echo "Run started at ${run_start_time}; Log directory will be ${Quaisar_H_log_directory}/${PROJECT}_on_${run_start_time}"
echo "Run started at ${run_start_time}" > "${log_dir}/${PROJECT}_on_${run_start_time}/${PROJECT}_on_${run_start_time}.log"
outarray+=("${PROJECT} started at ${run_start_time} and saved to ${PROJECT}_on_${run_start_time}.log")


#Each file in the list is put through the full pipeline
counter=1
for projfile in "${file_list[@]}";
do
	echo "${projfile}"
	file=$(echo "${projfile}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	proj=$(echo "${projfile}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	echo "${file} ${proj} ${BASEDIR}"
	if [[ -f "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" ]]; then
		rm "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"
	fi
	if [[ ! -f ${shareScript}/quaisar_${file}.sh ]]; then
		cp ${shareScript}/quaisar_template.sh ${shareScript}/quaisar_${file}.sh
		sed -i -e "s/quaisar_X/quaisar_${file}/g" "${shareScript}/quaisar_${file}.sh"
		sed -i -e "s/quasX/quasp_${file}/g" "${shareScript}/quaisar_${file}.sh"
		echo "Entering ${shareScript}/quaisar_${file}.sh" "${file}" "${proj}"
		qsub "${shareScript}/quaisar_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
		echo "Created and ran quaisar_${file}.sh"
	else
		echo "${shareScript}/quaisar_${file}.sh already exists, will resubmit"
		qsub "${shareScript}/quaisar_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
	fi
done

# Hold for completion of all submited single quaisars
for run_sample in "${file_list[@]}"; do
	waiting_sample=$(echo "${run_sample}" | cut -d'/' -f2)
	if [[ -f "${processed}/${PROJECT}/${waiting_sample}/${waiting_sample}_pipeline_stats.txt" ]]; then
		echo "${waiting_sample} is complete"
		mv ${shareScript}/quaisar_${waiting_sample}.sh ${log_dir}
		mv ${shareScript}/quaisar_${waiting_sample}.err ${log_dir}
		mv ${shareScript}/quaisar_${waiting_sample}.out ${log_dir}
	else
		while :
		do
				if [[ ${timer} -gt 1440 ]]; then
					echo "Timer exceeded limit of 86400 seconds = 24 hours"
					exit 1
				fi
				if [[ -f "${processed}/${PROJECT}/${waiting_sample}/${waiting_sample}_pipeline_stats.txt" ]]; then
					echo "${waiting_sample} is complete"
					mv ${shareScript}/quaisar_${waiting_sample}.sh ${log_dir}
					mv ${shareScript}/quaisar_${waiting_sample}.err ${log_dir}
					mv ${shareScript}/quaisar_${waiting_sample}.out ${log_dir}
					break
				else
					timer=$(( timer + 1 ))
					if [[ $(( timer % 5 )) -eq 0 ]]; then
						echo "Slept for ${timer} minutes so far"
					fi
					sleep 60
				fi
		done
	fi
done

# Concatenates lists if this run was an addition to an already processed folder
if [[ -f "${processed}/${PROJECT}/${PROJECT}_list_original.txt" ]]; then
	cat "${processed}/${PROJECT}/${PROJECT}_list.txt" >> "${processed}/${PROJECT}/${PROJECT}_list_original.txt"
	rm "${processed}/${PROJECT}/${PROJECT}_list.txt"
	mv "${processed}/${PROJECT}/${PROJECT}_list_original.txt" "${processed}/${PROJECT}/${PROJECT}_list.txt"
fi

runsumdate=$(date "+%m_%d_%Y_at_%Hh_%Mm")
${shareScript}/run_sum.sh ${PROJECT}
runsum=$(echo ${shareScript}/view_sum.sh ${PROJECT})
outarray+="${runsum}"

# Run the Seqlog creator on the proper file
if [ "${is_proj}" = "true" ]; then
	"${shareScript}/make_Seqlog_from_log.sh" "${PROJECT}"
else
	"${shareScript}/make_Seqlog_from_list.sh" "${processed}/${PROJECT}/${PROJECT}_list.txt"
fi

# Add print time the run completed in the text that will be emailed
global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Finished with run ${PROJECT} at ${global_end_time}"
outarray+=("
${PROJECT} finished at ${global_end_time}")
exit
#Send email to submitter and Nick with run status

if [ "${requestor}" != "nvx4" ]; then
	echo "Sending summary email to ${requestor}@cdc.gov & nvx4@cdc.gov"
	printf "%s\\n" "${outarray[@]}" | mail -s "Run Status for ${PROJECT}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
	printf "%s\\n" "${outarray[@]}" | mail -s "Run Status for ${PROJECT}_on_${run_start_time}_run.log" "${requestor}@cdc.gov"
else
	echo "Sending summary email to nvx4@cdc.gov"
	printf "%s\\n" "${outarray[@]}" | mail -s "Run Status for ${PROJECT}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
fi

# One final check for any dump files
if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
	echo "Found core dump files at end of processing ${PROJECT} and attempting to delete"
	find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
fi

if [[ -f "${shareScript}/config_${config_counter}.sh" ]]; then
	mv "${shareScript}/config_${config_counter}.sh" ${log_dir}
	:
fi

end_date=$(date "+%m_%d_%Y_at_%Hh_%Mm")
echo "Run ended at ${end_date}" >> "${log_dir}/${PROJECT}_on_${run_start_time}/${PROJECT}_on_${run_start_time}.log"
