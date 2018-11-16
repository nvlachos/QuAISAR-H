#!/bin/bash -l

#$ -o QuAISARX.out
#$ -e QuAISARX.err
#$ -N QuAISARX
#$ -cwd
#$ -q short.q

# Copy config file into working directory to allow changes to made to output directory if necessary
shareScript=$(pwd)
echo "${shareScript}"
if [[ -f config_template.sh ]]; then
	if [[ -f config.sh ]]; then
		rm -r config.sh
	fi
	echo "Trying to copy config_template.sh"
	cp config_template.sh config.sh
fi

#Import the config file with shortcuts and settings
. config.sh

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
. ${mod_changers}/list_modules.sh

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
project="${requestor}_${global_time}"
BASEDIR="${processed}"

for ((i=1 ; i <= nopts ; i++)); do
	#echo "${1} ${2}"
	case "${1}" in
		#Help/Usage section
		-h | --help)
			echo -e "\\n\\n\\n"
			echo "Usage: ./primary_processing.sh [i] [o] [d/n] [p/l/s] [project_name/list_of_samples.txt/sample_name]"
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
      #Gets output directory name of folder that all output files will be stored
      # NOT TESTED
      -o | --out-dir)
        BASEDIR="$2"
				project="$3"
        if [[ ! -d ${BASEDIR}/$3 ]]; then
            echo "Creating ${BASEDIR}/$3"
            mkdir -p ${BASEDIR}
        fi
        shift 3

				#processed=${BASEDIR}
        echo "processed=${BASEDIR}" >> "${shareScript}/config.sh"
				list_path="${BASEDIR}/${project}/${project}_list.txt"
        . ${shareScript}/config.sh
        echo "A - ${processed}"
        ;;
		#Gets name of folder that FASTA files will be in
		-i | --in-dir)
			INDATADIR="$2"
			if [[ -d  ${INDATADIR} ]]; then
				do_download="true"
				list_path="${BASEDIR}/${project}/${project}_list.txt"
			else
					echo "${INDATADIR} does not exist...exiting"
					exit 1
			fi
			postfix="$3"

			#echo "$INDATADIR $2"
			shift 3
			;;

		#Gets name of folder that FASTA files will be in
		-a | --assemblies-dir)
			INDATADIR="$2"
			if [[ -d  ${INDATADIR} ]]; then
				do_download="true"
				assembly_on="true"
				list_path="${BASEDIR}/${project}/${project}_list.txt"
			else
					echo "${INDATADIR} does not exist...exiting"
					exit 1
			fi
			postfix="$3"
			#echo "$INDATADIR $2"
			shift 3
			;;

		#Checks for (project) name of folder that all output files will be stored
		-p | --project-name)
			project="$2"
			if [[ -z ${INDATADIR} ]]; then
				INDATADIR="${processed}/${project}/"
			fi
			if [[ ${BASEDIR} = "${requestor}_${global_time}" ]]; then
				BASEDIR="${processed}"
			fi
			list_path="${processed}/${project}/${project}_list.txt"
			run_name="project_${project}"
			do_download="true"
			shift 2
			;;
		#Tells the script to run analyses from already downloaded fastq files
		-n | --no_download)
			do_download=false
			shift
			;;
		#Tells the script that only the files found in the attached list need to be run
		-lr | --list)
			list_path="$2"
			quick_list=$(echo "${2}" | cut -d'.' -f1)
			INDATADIR="${processed}/${project}" # NOT USED yet in list mode
			if [[ -z ${BASEDIR} ]]; then
				BASEDIR="${processed}"
			fi
			do_download="false"
			list_given="true"
			run_name="list_${quick_list}"
			shift 2
			;;
		#Tells the script that only the files found in the attached list need to be run
		-la | --list)
			list_path="$2"
			quick_list=$(echo "${2}" | cut -d'.' -f1)
			INDATADIR="${processed}/${project}" # NOT USED yet in list mode
			if [[ -z ${BASEDIR} ]]; then
				BASEDIR="${processed}"
			fi
			do_download="false"
			list_given="true"
			run_name="list_${quick_list}"
			shift 2
			;;
		#Captures any other characters in the args
		\?)
			echo "ERROR: ${BOLD}$2${NORM}is not a valid argument" >&2
			usage
			exit 1
			;;
	esac
done


# Short print out summary of run settings
echo -e "Source folder: ${INDATADIR}\\nOutput folder: ${BASEDIR}\\nDownload fastqs\(.gzs\): ${do_download}\\nList based analysis:  ${list_path}"








#project="nvx4_10-01-2018_at_12h_43m_08s"

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm")
#echo "${project} started at ${start_time} > ${processed}/${project}/${project}.log"
echo "${project} started at ${start_time}" > "${processed}/${project}/${project}.log"

do_assembly_download() {
	echo "Starting copying of Assemblies"
	check_time=$(date)
	#echo "Downloading started at ${check_time} >> ${processed}/${project}/${project}.log"
	echo "Downloading started at ${check_time}" >> "${processed}/${project}/${project}.log"
	mkdir /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}
	"${shareScript}/get_Assemblies_from_folder.sh" "${project}" "${INDATADIR}"
	echo "Done copying Assemblies from ${INDATADIR}"
	# Remove Undetermined zip files
	check_time=$(date)
	echo "Downloading finished at ${check_time}" >> "${processed}/${project}/${project}.log"
}


do_reads_download() {
	echo "Starting copying of FASTQs"
	check_time=$(date)
	#echo "Downloading started at ${check_time} >> ${processed}/${project}/${project}.log"
	echo "Downloading started at ${check_time}" >> "${processed}/${project}/${project}.log"
	mkdir ${processed}/${project}
	"${shareScript}/get_Reads_from_folder.sh" "${project}" "${INDATADIR}" "${postfix}"
	echo "Done copying FASTQs from input directory"
	# Remove Undetermined zip files
	#rm -r /scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/${project}/Undetermined*.gz
	check_time=$(date)
	echo "Downloading finished at ${check_time}" >> "${processed}/${project}/${project}.log"
}

declare sample_names
sample_index=0
main_dir="${processed}/${project}"
echo "${main_dir}"

#echo "> ${main_dir}/${project}_list.txt"
> "${main_dir}/${project}_list.txt"

#cp /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/config_template.sh ${main_dir}/config.sh
#echo "shareScript=/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts" >> ${main_dir}/config.sh



# Loop 1 - Make scripts to move and unzip all FASTQs to their individual isolate folders
make_fastq_unzipper() {
	for reads in ${BASEDIR}/${project}/*/FASTQs/*.gz; do
		echo "Loop 1 Make: ${reads}"

		if [[ "${reads}" = *"_R1_001"* ]]; then
			sample_ID_with_S=$(basename ${reads} _R1_001.fastq.gz)
			sample_ID=$(echo ${sample_ID_with_S} | rev | cut -d'_' -f2- | rev)
			echo -e "${reads}\n${sample_ID}"
			match=0
			for name in "${sample_names[@]}"; do
				if [[ "${name}" = "${sample_ID}" ]]; then
					echo "Name ${sample_ID} already exists in list (R1)"
					match=1
				else
					:
				fi
			done
			if [[ "${match}" = 0 ]]; then
				sample_names[${sample_index}]="${sample_ID}"
			fi

			# Create parent and FASTQ folders for current sample if they don't exist yet
			if [[ ! -d ${main_dir}/${sample_ID} ]]; then
				#echo "R1 - Creating ${main_dir}/${sample_ID}"
				mkdir ${main_dir}/${sample_ID}
			fi
			if [[ ! -d ${main_dir}/${sample_ID}/FASTQs ]]; then
				#echo "R1 - Creating ${main_dir}/${sample_ID}/FASTQs"
				mkdir ${main_dir}/${sample_ID}/FASTQs
			fi
			if [[ ! -d ${main_dir}/${sample_ID}/scripts ]]; then
				#echo "R1 - ${main_dir}/${sample_ID}/scripts"
				mkdir ${main_dir}/${sample_ID}/scripts
			fi
			if [[ ! -d ${main_dir}/${sample_ID}/logs ]]; then
				#echo "R1 - ${main_dir}/${sample_ID}/logs"
				mkdir ${main_dir}/${sample_ID}/logs
			fi

			#Create bash file for qsub submission to move and unzip isolate fastq.gz R1
			echo -e "#!/bin/bash -l\n" > "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "#$ -o getFASTQR1_${sample_ID}.out" >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "#$ -e getFASTQR1_${sample_ID}.err" >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "#$ -N getFASTQR1_${sample_ID}"   >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "echo $(date) > \"${main_dir}/${sample_ID}/logs/R1_unzipping_started.txt\"" >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "cp \"${reads}\" \"${main_dir}/${sample_ID}/FASTQs/${sample_ID}_R1_001.fastq.gz\"" >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "gunzip -c \"${main_dir}/${sample_ID}/FASTQs/${sample_ID}_R1_001.fastq.gz\" > \"${main_dir}/${sample_ID}/FASTQs/${sample_ID}_R1_001.fastq\"" >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "mv \"${reads}\" \"${main_dir}/${sample_ID}/FASTQs/${sample_ID}_R1_001.fastq.gz\"" >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"
			echo -e "echo $(date) > \"${main_dir}/${sample_ID}/logs/R1_unzipping_complete.txt\"" >> "${main_dir}/getFASTQR1_${sample_ID}_${start_time}.sh"


			# Add sample_ID to list of isolates to be run through other pipeline loops (backup to sample_names array)
			echo "${project}/${sample_ID}" >> "${main_dir}/${project}_list.txt"
			# Begin pipeline stats file for each sample by saving its original name
			full_name=$(basename ${reads} .gz)
			echo -e "---------Checking ${project}/${sample_ID} for successful completion----------" > "${main_dir}/${sample_ID}/${sample_ID}_pipeline_stats.txt"
			echo -e "Full identification: ${full_name}" >> "${main_dir}/${sample_ID}/${sample_ID}_pipeline_stats.txt"
			echo -e "Sample output folder: ${main_dir}/${project}/${sample_ID}/" >> "${main_dir}/${sample_ID}/${sample_ID}_pipeline_stats.txt"

			sample_index=$(( sample_index + 1 ))

		elif [[ "${reads}" = *"_R2_001"* ]]; then
			sample_ID_with_S=$(basename ${reads} _R2_001.fastq.gz)
			sample_ID=$(echo ${sample_ID_with_S} | rev | cut -d'_' -f2- | rev)
			match=0
			for name in "${sample_names[@]}"; do
				if [[ "${name}" = "${sample_ID}" ]]; then
					echo "Name ${sample_ID} already exists in list (R2)"
					match=1
				else
					:
				fi
			done
			if [[ "${match}" = 0 ]]; then
				sample_names[${sample_index}]="${sample_ID}"
			fi
			# Create parent and FASTQ folders for current sample if they don't exist yet
			if [[ ! -d ${main_dir}/${sample_ID} ]]; then
				#echo "R2 - Creating ${main_dir}/${sample_ID}"
				mkdir ${main_dir}/${sample_ID}
			fi
			if [[ ! -d ${main_dir}/${sample_ID}/FASTQs ]]; then
				#echo "R2 - Creating ${main_dir}/${sample_ID}/FASTQs"
				mkdir ${main_dir}/${sample_ID}/FASTQs
			fi
			if [[ ! -d ${main_dir}/${sample_ID}/scripts ]]; then
				#echo "R2 - Creating ${main_dir}/${sample_ID}/scripts"
				mkdir ${main_dir}/${sample_ID}/scripts
			fi
			if [[ ! -d ${main_dir}/${sample_ID}/logs ]]; then
				#echo "R1 - ${main_dir}/${sample_ID}/logs"
				mkdir ${main_dir}/${sample_ID}/logs
			fi
			#Create bash file for qsub submission to move and unzip isolate fastq.gz R1
			echo -e "#!/bin/bash -l\n" > "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "#$ -o getFASTQR2_${sample_ID}.out" >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "#$ -e getFASTQR2_${sample_ID}.err" >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "#$ -N getFASTQR2_${sample_ID}"   >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "echo $(date) > \"${main_dir}/${sample_ID}/logs/R2_unzipping_started.txt\"" >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "cp \"${reads}\" \"${main_dir}/${sample_ID}/FASTQs/${sample_ID}_R2_001.fastq.gz\"" >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "gunzip -c \"${main_dir}/${sample_ID}/FASTQs/${sample_ID}_R2_001.fastq.gz\" > \"${main_dir}/${sample_ID}/FASTQs/${sample_ID}_R2_001.fastq\"" >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "mv \"${reads}\" \"${main_dir}/${sample_ID}/FASTQs/${sample_ID}_R2_001.fastq.gz\"" >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
			echo -e "echo $(date) > \"${main_dir}/${sample_ID}/logs/R2_unzipping_complete.txt\"" >> "${main_dir}/getFASTQR2_${sample_ID}_${start_time}.sh"
		else
			echo "Unknown .gz file: ${reads}"
		fi
	done
}

# Makes a list of samples from folders available in
make_list_from_folder() {
	#echo "Trying - ${main_dir}/"
	counter=0
	for dir in ${main_dir}/*/
	do
		isolate=$(basename ${dir})
		#echo "${isolate}"
		sample_names[${counter}]="${isolate}"
		echo "${project}/${isolate}" >> "${main_dir}/${project}_list.txt"
		counter=$(( counter + 1 ))
		#if [[ -d "${dir}" ]]; then
		#	echo "${dir}"
		#fi
	done

	echo "${sample_names[@]}"
}

make_list_from_list() {
	counter=0
	while IFS= read -r var; do
		project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
		sample_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
		sample_names[${counter}]=${sample_name}
	done < ${list_path}
}

# # Submit Loop 1 scripts
submit_fastq_unzipper() {
	check_time=$(date)
	echo "Loop 1 started at ${check_time}" >> "${processed}/${project}/${project}.log"
	for sample in "${sample_names[@]}"; do
		echo "Loop 1 Send: ${sample}"

		echo -e 'Trying to call L1-Unzip/Move qsubs' #\n${main_dir}/getFASTQR1_${sample}_${start_time}.sh\n${main_dir}/getFASTQR2_${sample}_${start_time}.sh'
		qsub "${main_dir}/getFASTQR1_${sample}_${start_time}.sh"
		qsub "${main_dir}/getFASTQR2_${sample}_${start_time}.sh"
	done
	check_time=$(date)
	echo "Loop 1 end at ${check_time}" >> "${processed}/${project}/${project}.log"
}

# Loop 2 Make script files for all tools requiring unzipped FASTQ files
make_relies_on_unzipped_fastqs() {
	for sample in "${sample_names[@]}"
	do
		echo "Loop 2 Make: ${sample}"
		# Check for proper pairing later

		#Create script to run QC on reads
		echo -e "#!/bin/bash -l\n" > "${main_dir}/QC_${sample}_${start_time}.sh"
		echo -e "#$ -o QC_${sample}.out" >> "${main_dir}/QC_${sample}_${start_time}.sh"
		echo -e "#$ -e QC_${sample}.err" >> "${main_dir}/QC_${sample}_${start_time}.sh"
		echo -e "#$ -N QC_${sample}"   >> "${main_dir}/QC_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/QC_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/QC_${sample}_${start_time}.sh"
		echo -e "module load Python/2.7.15\n" >> "${main_dir}/QC_${sample}_${start_time}.sh"
		echo -e "mkdir \"${main_dir}/${sample}/preQCcounts\"" >> "${main_dir}/QC_${sample}_${start_time}.sh"
		echo -e "python \"${shareScript}/Fastq_Quality_Printer.py\" \"${main_dir}/${sample}/FASTQs/${sample}_R1_001.fastq\" \"${main_dir}/${sample}/FASTQs/${sample}_R2_001.fastq\" > \"${main_dir}/${sample}/preQCcounts/${sample}_counts.txt\"" >> "${main_dir}/QC_${sample}_${start_time}.sh"

		#Create script to run BBDuk, Trimmomatic, and QC on trimmed reads
		echo -e "#!/bin/bash -l\n" > "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "#$ -o BTQC_${sample}.out" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "#$ -e BTQC_${sample}.err" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "#$ -N BTQC_${sample}"   >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "#$ -pe smp ${procs}" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "module load BBMap/35.92" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "module load trimmomatic/0.36" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "module load Python/2.7.15" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_trimming_started.txt\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "mkdir \"${main_dir}/${sample}/removedAdapters\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "mkdir \"${main_dir}/${sample}/trimmed\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "bbduk.sh -\"${bbduk_mem}\" threads=\"${procs}\" in=\"${main_dir}/${sample}/FASTQs/${sample}_R1_001.fastq\" in2=\"${main_dir}/${sample}/FASTQs/${sample}_R2_001.fastq\" out=\"${main_dir}/${sample}/removedAdapters/${sample}-noPhiX-R1.fsq\" out2=\"${main_dir}/${sample}/removedAdapters/${sample}-noPhiX-R2.fsq\" ref=\"${local_DBs}/phiX.fasta\" k=\"${bbduk_k}\" hdist=\"${bbduk_hdist}\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "trimmomatic \"${trim_endtype}\" -\"${trim_phred}\" -threads \"${procs}\" \"${main_dir}/${sample}/removedAdapters/${sample}-noPhiX-R1.fsq\" \"${main_dir}/${sample}/removedAdapters/${sample}-noPhiX-R2.fsq\" \"${main_dir}/${sample}/trimmed/${sample}_R1_001.paired.fq\" \"${main_dir}/${sample}/trimmed/${sample}_R1_001.unpaired.fq\" \"${main_dir}/${sample}/trimmed/${sample}_R2_001.paired.fq\" \"${main_dir}/${sample}/trimmed/${sample}_R2_001.unpaired.fq\" ILLUMINACLIP:\"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}\" SLIDINGWINDOW:\"${trim_window_size}:${trim_window_qual}\" LEADING:\"${trim_leading}\" TRAILING:\"${trim_trailing}\" MINLEN:\"${trim_min_length}\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "python \"${shareScript}/Fastq_Quality_Printer.py\" \"${main_dir}/${sample}/trimmed/${sample}_R1_001.paired.fq\" \"${main_dir}/${sample}/trimmed/${sample}_R2_001.paired.fq\" > \"${main_dir}/${sample}/preQCcounts/${sample}_trimmed_counts.txt\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		# Merge both unpaired fq files into one for GOTTCHA
		echo -e "cat \"${main_dir}/${sample}/trimmed/${sample}_R1_001.paired.fq\" \"${main_dir}/${sample}/trimmed/${sample}_R2_001.paired.fq\" > \"${main_dir}/${sample}/trimmed/${sample}.paired.fq\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "cat \"${main_dir}/${sample}/trimmed/${sample}_R1_001.unpaired.fq\" \"${main_dir}/${sample}/trimmed/${sample}_R2_001.unpaired.fq\" > \"${main_dir}/${sample}/trimmed/${sample}.single.fq\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "gzip < \"${main_dir}/${sample}/trimmed/${sample}_R1_001.paired.fq\" > \"${main_dir}/${sample}/trimmed/${sample}_R1_001.paired.fq.gz\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "gzip < \"${main_dir}/${sample}/trimmed/${sample}_R2_001.paired.fq\" > \"${main_dir}/${sample}/trimmed/${sample}_R2_001.paired.fq.gz\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_trimming_complete.txt\"" >> "${main_dir}/BTQC_${sample}_${start_time}.sh"
	done
}

# Loop 2 - Submit all scripts for tools requiring unzipped FASTQs, ensuring FASTQs are available to work with
submit_relies_on_unzipped_fastqs() {
	check_time=$(date)
	echo "Loop 2 started at ${check_time}" >> "${processed}/${project}/${project}.log"
	for sample in "${sample_names[@]}";
	do
		echo "Loop 2 Send: ${sample}"
		unzip_success="false"
		counter=0
		total=0
		while :
		do
			if [[ -f ${main_dir}/${sample}/logs/R1_unzipping_complete.txt ]] && [[ -f ${main_dir}/${sample}/logs/R2_unzipping_complete.txt ]]; then
				unzip_success="true"
				break
			elif [[ ${counter} -gt 360 ]]; then
				break
			else
				counter=$(( counter + 1 ))
				total=$(( counter * 5 ))
				echo "Waiting 5 seconds (total=${total}s) for ${main_dir}/${sample}/logs/[R1/R2]]_ziping_complete(s).txt to appear"
				sleep 5s
			fi
		done
		if [[ "${unzip_success}" = "true" ]]; then
			echo -e 'Cleaning L1-Unzip/Move qsubs'
			mv "${main_dir}/getFASTQR1_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/"
			mv "${shareScript}/getFASTQR1_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/getFASTQR2_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/"
			mv "${shareScript}/getFASTQR2_${sample}"* "${main_dir}/${sample}/scripts/"
			# Trying to submit QC and Trimming scripts"
			echo -e 'Trying to call L2-QCBT qsubs' #\n${main_dir}/QC_${sample}_${start_time}.sh\n${main_dir}/BTQC_${sample}_${start_time}.sh'
			qsub "${main_dir}/QC_${sample}_${start_time}.sh"
			qsub "${main_dir}/BTQC_${sample}_${start_time}.sh"
		else
			echo "30 minutes has elapsed and ${main_dir}/${sample}/FASTQs/${sample}_L001_R1_001.fastq has not appeared"
		fi
	done
	check_time=$(date)
	echo "Loop 2 ended at ${check_time}" >> "${processed}/${project}/${project}.log"
}

# Loop 3 - Make scripts for all tools requiring trimmed reads
make_relies_on_trimmed_fastqs() {
	for sample in "${sample_names[@]}"
	do
		echo "Loop 3 Make: ${sample}"

		# # Create script to run Kraken on raw reads
		echo -e "#!/bin/bash -l\n" > "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "#$ -o krakr_${sample}.out" >> "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "#$ -e krakr_${sample}.err" >> "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "#$ -N krakr_${sample}"   >> "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "#$ -pe smp ${procs}" >> "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "module load kraken/1.0.0" >> "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "module load krona/2.7\n" >> "${main_dir}/krakr_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_kraken.sh\" \"${sample}\" pre paired \"${project}\"" >> "${main_dir}/krakr_${sample}_${start_time}.sh"

		# Create script to run GOTTCHA
		echo -e "#!/bin/bash -l\n" > "${main_dir}/gott_${sample}_${start_time}.sh"
		echo -e "#$ -o gott_${sample}.out" >> "${main_dir}/gott_${sample}_${start_time}.sh"
		echo -e "#$ -e gott_${sample}.err" >> "${main_dir}/gott_${sample}_${start_time}.sh"
		echo -e "#$ -N gott_${sample}"   >> "${main_dir}/gott_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/gott_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/gott_${sample}_${start_time}.sh"
		echo -e "module load gottcha\n" >> "${main_dir}/gott_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_gottcha.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/gott_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_GOTT_complete.txt\"" >> "${main_dir}/gott_${sample}_${start_time}.sh"

		# Create script to call SRST2 AR
		echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		echo -e "module load srst2/0.1.7" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		echo -e "module unload Python/2.7.11" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		echo -e "module load Python/2.7.15\n" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
		echo -e "\"${shareScript}/run_srst2_on_singleDB.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"

		# Create script to call normal SPAdes
		echo -e "#!/bin/bash -l\n" > "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "#$ -o SPAdn_${sample}.out" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "#$ -e SPAdn_${sample}.err" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "#$ -N SPAdn_${sample}"   >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "#$ -pe smp ${procs}" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "module load SPAdes/3.12.0\n" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_SPAdes.sh\" \"${sample}\" normal \"${project}\"" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_full_assembling_complete.txt\"" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "python \"${shareScript}/removeShortContigs.py\" \"${main_dir}/${sample}/Assembly/scaffolds.fasta\" 500" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "mv \"${main_dir}/${sample}/Assembly/scaffolds.fasta.TRIMMED.fasta\" \"${main_dir}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta\"" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_full_assembly_trimming_complete.txt\"" >> "${main_dir}/SPAdn_${sample}_${start_time}.sh"

		# Create script to call plasmid SPAdes
		echo -e "#!/bin/bash -l\n" > "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "#$ -o SPAdp_${sample}.out" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "#$ -e SPAdp_${sample}.err" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "#$ -N SPAdp_${sample}"   >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "#$ -pe smp ${procs}" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "module load SPAdes/3.12.0\n" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_SPAdes.sh\" \"${sample}\" plasmid \"${project}\"" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_plasmid_assembling_complete.txt\"" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "python \"${shareScript}/removeShortContigs.py\" \"${main_dir}/${sample}/plasmidAssembly/scaffolds.fasta\" 500" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "mv \"${main_dir}/${sample}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta\" \"${main_dir}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta\"" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_plasmid_assembly_trimming_complete.txt\"" >> "${main_dir}/SPAdp_${sample}_${start_time}.sh"
	done
}

#Loop 3 - Submit all scripts for Tools requiring trimmed fastqs
submit_relies_on_trimmed_fastqs() {
	check_time=$(date)
	echo "Loop 3 started at ${check_time}" >> "${processed}/${project}/${project}.log"
	for sample in "${sample_names[@]}";
	do
		echo "Loop 3 Send: ${sample}"
		trim_success="false"
		counter=0
		total=0
		while :
		do
			if [[ -f "${main_dir}/${sample}/logs/${sample}_trimming_complete.txt" ]]; then
				trim_success="true"
				break
			elif [[ ${counter} -gt 360 ]]; then
				break
			else
				counter=$(( counter + 1 ))
				total=$(( counter * 5 ))
				echo "Waiting 5 seconds (total=${total}s) for ${main_dir}/${sample}/logs/${sample}_trimming_complete.txt to appear"
				sleep 5
			fi
		done
		if [[ "${trim_success}" = "true" ]]; then
			:
			# Cleanup QC and BTQC scripts
			echo 'Cleaning L2-QCBT scripts'
			mv "${main_dir}/BTQC_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/"
			mv "${shareScript}/BTQC_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/QC_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/"
			mv "${shareScript}/QC_${sample}"* "${main_dir}/${sample}/scripts/"
			echo -e 'Trying to call L3-paired.fq_reliant qsubs' #\n${main_dir}/krakr_${sample}_${start_time}.sh\n${main_dir}/gott_${sample}_${start_time}.sh\n${main_dir}/srst2AR_${sample}_${start_time}.sh\n${main_dir}/SPAdn_${sample}_${start_time}.sh\n${main_dir}/SPAdp_${sample}_${start_time}.sh'
			qsub "${main_dir}/krakr_${sample}_${start_time}.sh"
			qsub "${main_dir}/gott_${sample}_${start_time}.sh"
			qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			qsub "${main_dir}/SPAdn_${sample}_${start_time}.sh"
			qsub "${main_dir}/SPAdp_${sample}_${start_time}.sh"
		else
			echo "30 minutes has elapsed and ${main_dir}/${sample}/trimmed/${sample}_R1_001.paired.fq or ${main_dir}/${sample}/trimmed/${sample}_R2_001.paired.fq has not appeared"
		fi
	done
	check_time=$(date)
	echo "Loop 3 ended at ${check_time}" >> "${processed}/${project}/${project}.log"
}

# Loop 3.5 - Check that SPAdes completed successfully
check_for_Assemblies() {
	for sample in "${sample_name[@]}";
	do
		if [[ ! -s "${main_dir}/${sample}/Assembly/scaffolds.fasta" ]]; then
			# Create script to call normal SPAdes
			#re qusb Assembly
			echo "$sample did not finish normal assembling"
			if [[ -f "${main_dir}/${sample}/logs/${sample}_full_assembling_complete.txt" ]]; then
				echo "Says it finished in logs"
			else
				echo "Says it did not finish in logs"
			fi
		else
			echo "${sample} scaffolds does not exist or is 0 size"
		fi
		if [[ ! -s "${main_dir}/${sample}/plasmidAssembly/scaffolds.fasta" ]]; then
			# Create script to call plasmid SPAdes
			#re qsub plasmidAssembly
			echo "$sample did not finish plamsid assembling"
			if [[ -f "${main_dir}/${sample}/logs/${sample}_plasmid_assembing_complete.txt" ]]; then
				echo "Says it finished in logs"
			else
				echo "Says it did not finish in logs"
			fi
		else
			echo "${sample} plasmid_scaffolds does not exist or is 0 size"
		fi
	done
}

# Loop 4 - Make scripts for tools requiring trimmed assembly
make_relies_on_trimmed_assemblies() {
	for sample in "${sample_names[@]}";
	do
		echo "Loop 4 Make: ${sample}"

		# Create script to call QUAST
		echo -e "#!/bin/bash -l\n" > "${main_dir}/QUAST_${sample}_${start_time}.sh"
		echo -e "#$ -o QUAST_${sample}.out" >> "${main_dir}/QUAST_${sample}_${start_time}.sh"
		echo -e "#$ -e QUAST_${sample}.err" >> "${main_dir}/QUAST_${sample}_${start_time}.sh"
		echo -e "#$ -N QUAST_${sample}"   >> "${main_dir}/QUAST_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/QUAST_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/QUAST_${sample}_${start_time}.sh"
		echo -e "module load quast/4.3" >> "${main_dir}/QUAST_${sample}_${start_time}.sh"
		echo -e "module load Python/2.7.15" >> "${main_dir}/QUAST_${sample}_${start_time}.sh"
		# No modules needed to run QUAST
		echo -e "\"${shareScript}/run_Assembly_Quality_Check.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/QUAST_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_QUAST_complete.txt\"" >> "${main_dir}/QUAST_${sample}_${start_time}.sh"

		# Create scripts for Prokka
		echo -e "#!/bin/bash -l\n" > "${main_dir}/PROKK_${sample}_${start_time}.sh"
		echo -e "#$ -o PROKK_${sample}.out" >> "${main_dir}/PROKK_${sample}_${start_time}.sh"
		echo -e "#$ -e PROKK_${sample}.err" >> "${main_dir}/PROKK_${sample}_${start_time}.sh"
		echo -e "#$ -N PROKK_${sample}"   >> "${main_dir}/PROKK_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/PROKK_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/PROKK_${sample}_${start_time}.sh"
		echo -e "module load prokka/1.12\n" >> "${main_dir}/PROKK_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_prokka.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/PROKK_${sample}_${start_time}.sh"
		echo -e "mv \"${main_dir}/${sample}/Assembly/scaffolds.fasta.TRIMMED.fasta\" \"${main_dir}/${sample}/Assembly/scaffolds.fasta.TRIMMED_original.fasta\""
		echo -e "\"python3 ${shareScript}/fasta_headers.py\" \"${main_dir}/${sample}/Assembly/scaffolds.fasta.TRIMMED_original.fasta\" \"${main_dir}/${sample}/Assembly/scaffolds.fasta.TRIMMED.fasta\""
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_PROKK_complete.txt\"" >> "${main_dir}/PROKK_${sample}_${start_time}.sh"

		# Create scripts for c-SSTAR
		echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
		echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
		echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
		echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
		echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
		# Defaulting to gapped/98, change if you want to include user preferences
		echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g h \"${project}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"

		# Create script to run plasmidFinder on assembly
		echo -e "#!/bin/bash -l\n" > "${main_dir}/pFinf_${sample}_${start_time}.sh"
		echo -e "#$ -o pFinf_${sample}.out" >> "${main_dir}/pFinf_${sample}_${start_time}.sh"
		echo -e "#$ -e pFinf_${sample}.err" >> "${main_dir}/pFinf_${sample}_${start_time}.sh"
		echo -e "#$ -N pFinf_${sample}"   >> "${main_dir}/pFinf_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/pFinf_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/pFinf_${sample}_${start_time}.sh"
		echo -e "module load PlasmidFinder/1.3\n" >> "${main_dir}/pFinf_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_plasmidFinder.sh\" \"${sample}\" \"${project}\" \"plasmid\"" >> "${main_dir}/pFinf_${sample}_${start_time}.sh"

		# Create script to run MLST on Assembly
		echo -e "#!/bin/bash -l\n" > "${main_dir}/MLST_${sample}_${start_time}.sh"
		echo -e "#$ -o MLST_${sample}.out" >> "${main_dir}/MLST_${sample}_${start_time}.sh"
		echo -e "#$ -e MLST_${sample}.err" >> "${main_dir}/MLST_${sample}_${start_time}.sh"
		echo -e "#$ -N MLST_${sample}"   >> "${main_dir}/MLST_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/MLST_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/MLST_${sample}_${start_time}.sh"
		echo -e "module load mlst/2.9\n" >> "${main_dir}/MLST_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_MLST.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/MLST_${sample}_${start_time}.sh"

		# Create script to run blast16s ID on assembly
		echo -e "#!/bin/bash -l\n" > "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "#$ -o blast16sID_${sample}.out" >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "#$ -e blast16sID_${sample}.err" >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "#$ -N blast16sID_${sample}"   >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "#$ -pe smp ${procs}" >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "module load barrnap/0.8\n" >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/16s_blast.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_blast16sID_complete.txt\"" >> "${main_dir}/blast16sID_${sample}_${start_time}.sh"

		# Create script to run Kraken on assembly
		echo -e "#!/bin/bash -l\n" > "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "#$ -o kraka_${sample}.out" >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "#$ -e kraka_${sample}.err" >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "#$ -N kraka_${sample}"   >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "#$ -pe smp ${procs}" >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "module load kraken/1.0.0" >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "module load krona/2.7" >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "module load Python/2.7.15\n" >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_kraken.sh\" \"${sample}\" post assembled \"${project}\"" >> "${main_dir}/kraka_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_kraka_complete.txt\"" >> "${main_dir}/kraka_${sample}_${start_time}.sh"

		# Create script to run Kraken on assembly (FULL Database)
		echo -e "#!/bin/bash -l\n" > "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "#$ -o krakf_${sample}.out" >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "#$ -e krakf_${sample}.err" >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "#$ -N krakf_${sample}"   >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "#$ -pe smp ${procs}" >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "module load kraken/1.0.0" >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "module load krona/2.7" >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "module load Python/2.7.15\n" >> "${main_dir}/krakf_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_kraken.sh\" \"${sample}\" post assembled \"${project}\" full" >> "${main_dir}/krakf_${sample}_${start_time}.sh"

		# Make script to run ANI
		echo -e "#!/bin/bash -l\n" > "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "#$ -o ANI_${sample}.out" >> "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "#$ -e ANI_${sample}.err" >> "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "#$ -N ANI_${sample}"   >> "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "#$ -pe smp ${procs}" >> "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "module load pyani\n" >> "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_ANI.sh\" \"${sample}\" All All \"${project}\"" >> "${main_dir}/ANI_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_ANI_complete.txt\"" >> "${main_dir}/ANI_${sample}_${start_time}.sh"

		# Make script to run PlasFlow
		echo -e "#!/bin/bash -l\n" > "${main_dir}/plasFlow_${sample}_${start_time}.sh"
		echo -e "#$ -o plasFlow_${sample}.out" >> "${main_dir}/plasFlow_${sample}_${start_time}.sh"
		echo -e "#$ -e plasFlow_${sample}.err" >> "${main_dir}/plasFlow_${sample}_${start_time}.sh"
		echo -e "#$ -N plasFlow_${sample}"   >> "${main_dir}/plasFlow_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/plasFlow_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/plasFlow_${sample}_${start_time}.sh"
		echo -e "module load PlasFlow/1.0\n" >> "${main_dir}/plasFlow_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_plasFlow.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/plasFlow_${sample}_${start_time}.sh"

		# Check if plasmid SPAdes produced any output...if so, make scripts to run csstar and plasmidFinder on it
		if [[ -s "${main_dir}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta}" ]]; then
			echo -e "mv \"${main_dir}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds.fasta.TRIMMED.fasta\" \"${main_dir}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds.fasta.TRIMMED_original.fasta\""
			echo -e "\"python3 ${shareScript}/fasta_headers.py\" \"${main_dir}/${sample}/Assembly/${sample}_plasmid_scaffolds.fasta.TRIMMED_original.fasta\" \"${main_dir}/${sample}/Assembly/${sample}_plasmid_scaffolds.fasta.TRIMMED.fasta\""

			echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
			echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
			echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
			echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
			echo - "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
			# Defaulting to gapped/40, change if you want to include user preferences
			echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g o \"${project}\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"

			echo -e "#!/bin/bash -l\n" > "${main_dir}/pFinp_${sample}_${start_time}.sh"
			echo -e "#$ -o pFinp_${sample}.out" >> "${main_dir}/pFinp_${sample}_${start_time}.sh"
			echo -e "#$ -e pFinp_${sample}.err" >> "${main_dir}/pFinp_${sample}_${start_time}.sh"
			echo -e "#$ -N pFinp_${sample}"   >> "${main_dir}/pFinp_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/pFinp_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/pFinp_${sample}_${start_time}.sh"
			echo -e "module load PlasmidFinder/1.3\n" >> "${main_dir}/pFinp_${sample}_${start_time}.sh"
			echo -e "\"${shareScript}/run_plasmidFinder.sh\" \"${sample}\" \"${project}\" \"plasmid_on_plasmidAssembly\"" >> "${main_dir}/pFinp_${sample}_${start_time}.sh"
		fi
	done
}

# Loop 4 - Submit scripts for tools requiring trimmed assembly
submit_relies_on_trimmed_assemblies() {
	check_time=$(date)
	echo "Loop 4 started at ${check_time}" >> "${processed}/${project}/${project}.log"
	for sample in "${sample_names[@]}";
	do
		echo "Loop 4 Send: ${sample}"
		assembly_success="false"
		counter=0
		total=0
		while :
		do
			if [[ -f "${main_dir}/${sample}/logs/${sample}_plasmid_assembly_trimming_complete.txt" ]] && [[ -f "${main_dir}/${sample}/logs/${sample}_full_assembly_trimming_complete.txt" ]]; then
				assembly_success="true"
				break
			elif [[ "${counter}" -gt 360 ]]; then
				break
			else
				counter=$(( counter + 1 ))
				total=$(( counter * 5 ))
				echo "Waiting 5 seconds (total=${total}s) for ${main_dir}/${sample}/logs/${sample}_(full/&/plasmid)_assembly_trimming_complete.txt to appear"
				sleep 5
			fi
		done
		if [[ "${assembly_success}" = "true" ]]; then
			:
			# Cleanup L3 (paired.fq reliant) scripts
			echo 'Cleaning L3-paired.fq_reliant scripts'
			mv "${main_dir}/krakr_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/krakr_${sample}_${start_time}.sh"
			mv "${shareScript}/krakr_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/gott_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/gott_${sample}_${start_time}.sh"
			mv "${shareScript}/gott_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/srst2AR_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/srst2AR_${sample}_${start_time}.sh"
			mv "${shareScript}/srst2AR_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/SPAdn_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/SPAdn_${sample}_${start_time}.sh"
			mv "${shareScript}/SPAdn_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/SPAdp_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/SPAdp_${sample}_${start_time}.sh"
			mv "${shareScript}/SPAdp_${sample}"* "${main_dir}/${sample}/scripts/"
			echo -e 'Trying to call L4-trimmed_assembly_reliant qsubs' #\n${main_dir}/QUAST_${sample}_${start_time}.sh\n${main_dir}/kraka_${sample}_${start_time}.sh\n${main_dir}/PROKK_${sample}_${start_time}.sh\n${main_dir}/MLST_${sample}_${start_time}.sh\n${main_dir}/blast16sID_${sample}_${start_time}.sh\n${main_dir}/csstn_${sample}_${start_time}.sh\n${main_dir}/pFinf_${sample}_${start_time}.sh\n${main_dir}/ANI_${sample}_${start_time}.sh\n${main_dir}/plasFlow_${sample}_${start_time}.sh\n(O)${main_dir}/pFinp_${sample}_${start_time}.sh\n(O)${main_dir}/csstp_${sample}_${start_time}.sh'
			qsub "${main_dir}/QUAST_${sample}_${start_time}.sh"
			# Should we put a check here to prevent unnecessary calls for bad samples???
			qsub "${main_dir}/kraka_${sample}_${start_time}.sh"
			qsub "${main_dir}/PROKK_${sample}_${start_time}.sh"
			qsub "${main_dir}/MLST_${sample}_${start_time}.sh"
			qsub "${main_dir}/blast16sID_${sample}_${start_time}.sh"
			qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
			qsub "${main_dir}/pFinf_${sample}_${start_time}.sh"
			qsub "${main_dir}/ANI_${sample}_${start_time}.sh"
#			qsub "${main_dir}/plasFlow_${sample}_${start_time}.sh"
			if [[ -s "${main_dir}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta}" ]]; then
				qsub "${main_dir}/pFinp_${sample}_${start_time}.sh"
				qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
			fi
		else
			echo "30 minutes has elapsed and ${main_dir}/${sample}/logs/${sample}_plasmid_assembly_trimming_complete.txt or ${main_dir}/${sample}/logs/${sample}_full_assembly_trimming_complete.txt has not appeared"
		fi
	done
	check_time=$(date)
	echo "Loop 4 ended at ${check_time}" >> "${processed}/${project}/${project}.log"
}

# Loop 5 - Make script to provide confirmation of species for each sample
make_relies_on_species_files() {
	for sample in "${sample_names[@]}";
	do
		echo "Loop 5 Make: ${sample}"
		echo -e "#!/bin/bash -l\n" > "${main_dir}/taxID_${sample}_${start_time}.sh"
		echo -e "#$ -o taxID_${sample}.out" >> "${main_dir}/taxID_${sample}_${start_time}.sh"
		echo -e "#$ -e taxID_${sample}.err" >> "${main_dir}/taxID_${sample}_${start_time}.sh"
		echo -e "#$ -N taxID_${sample}"   >> "${main_dir}/taxID_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/taxID_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/taxID_${sample}_${start_time}.sh"
		# No modules needed to run determine_taxID
		echo -e "\"${shareScript}/determine_taxID.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/taxID_${sample}_${start_time}.sh"
	done
}

# Loop 5 - Submit script confirming species identification
submit_relies_on_species_files() {
	check_time=$(date)
	echo "Loop 5 started at ${check_time}" >> "${processed}/${project}/${project}.log"
	for sample in "${sample_names[@]}";
	do
		echo "Loop 5 Send: ${sample}"
		counter=0
		total=0
		species_success="false"
		while :
		do
			#if [[ -f "${processed}/${project}/${sample}/${sample}.tax" ]]; then
			if [[ -f "${main_dir}/${sample}/logs/${sample}_ANI_complete.txt" ]] && [[ -f "${main_dir}/${sample}/logs/${sample}_kraka_complete.txt" ]] && [[ -f "${main_dir}/${sample}/logs/${sample}_GOTT_complete.txt" ]] && [[ -f "${main_dir}/${sample}/logs/${sample}_blast16sID_complete.txt" ]]; then
				 species_success="true"
				 break
			elif [[ ${counter} -gt 360 ]]; then
				if [[ -f "${main_dir}/${sample}/logs/${sample}_ANI_complete.txt" ]] || [[ -f "${main_dir}/${sample}/logs/${sample}_kraka_complete.txt" ]] || [[ -f "${main_dir}/${sample}/logs/${sample}_GOTT_complete.txt" ]] || [[ -f "${main_dir}/${sample}/logs/${sample}_blast16sID_complete.txt" ]]; then
					species_success="true"
					break
				else
					break
				fi
			else
				counter=$(( counter + 1 ))
				total=$(( counter * 5 ))
				echo "Waiting 5 seconds (total=${total}s) for ALL SPECIES files to appear"
				sleep 5
			fi
		done
		if [[ "${species_success}" == "true" ]]; then
			#Clean Loop 4 script files
			echo 'Cleaning L4-assembly_reliant scripts'
			mv "${main_dir}/QUAST_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/QUAST_${sample}_${start_time}.sh"
			mv "${shareScript}/QUAST_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/kraka_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/kraka_${sample}_${start_time}.sh"
			mv "${shareScript}/kraka_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/PROKK_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/PROKK_${sample}_${start_time}.sh"
			mv "${shareScript}/PROKK_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/krakf_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/krakf_${sample}_${start_time}.sh"
			mv "${shareScript}/krakf_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/MLST_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/MLST_${sample}_${start_time}.sh"
			mv "${shareScript}/MLST_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/blast16sID_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/blast16sID_${sample}_${start_time}.sh"
			mv "${shareScript}/blast16sID_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/csstn_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/ccstn_${sample}_${start_time}.sh"
			mv "${shareScript}/csstn_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/pFinf_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/pFinf_${sample}_${start_time}.sh"
			mv "${shareScript}/pFinf_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/ANI_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/ANI_${sample}_${start_time}.sh"
			mv "${shareScript}/ANI_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/plasFlow_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/plasFlow_${sample}_${start_time}.sh"
			mv "${shareScript}/plasFlow_${sample}"* "${main_dir}/${sample}/scripts/"
			if [[ -s "${main_dir}/${sample}/plasmidAssembly/${sample}_plasmid_scaffolds_trimmed.fasta}" ]]; then
				mv "${main_dir}/pFinp_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/pFinp_${sample}_${start_time}.sh"
				mv "${shareScript}/pFinp_${sample}"* "${main_dir}/${sample}/scripts/"
				mv "${main_dir}/csstp_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/csstp_${sample}_${start_time}.sh"
				mv "${shareScript}/csstp_${sample}"* "${main_dir}/${sample}/scripts/"
			fi
			echo -e 'Trying to call L5-species_determining qsubs' #\n${main_dir}/taxID_${sample}_${start_time}.sh'
			qsub "${main_dir}/taxID_${sample}_${start_time}.sh"
		else
			echo "None of the taxonomy files appeared within 30 minutes..."
		fi
	done
	check_time=$(date)
	echo "Loop 5 ended at ${check_time}" >> "${processed}/${project}/${project}.log"
}

# Loop 6 - Make scripts for tools requiring species identification
make_relies_on_species_confirmation() {
	for sample in "${sample_names[@]}";
	do
		echo "Loop 6 Make: ${sample}"
		species_success="false"
		counter=0
		total=0
		while :
		do
			if [[ -f "${main_dir}/${sample}/${sample}.tax" ]]; then
				species_success="true"
				while IFS= read -r line;
				do
						# Grab first letter of line (indicating taxonomic level)
						first=${line::1}
						echo ${line}
						# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
						if [ "${first}" = "s" ]
						then
							species=$(echo "${line}" | cut -d'	' -f2 | tr -d [:space:])
						elif [ "${first}" = "G" ]
						then
							genus=$(echo "${line}" | cut -d'	' -f2 | tr -d [:space:])
						elif [ "${first}" = "F" ]
						then
							family=$(echo "${line}" | cut -d'	' -f2 | tr -d [:space:])
						elif [ "${first}" = "O" ]
						then
							order=$(echo "${line}" | cut -d'	' -f2 | tr -d [:space:])
						elif [ "${first}" = "C" ]
						then
							class=$(echo "${line}" | cut -d'	' -f2 | tr -d [:space:])
						elif [ "${first}" = "P" ]
						then
							phylum=$(echo "${line}" | cut -d'	' -f2 | tr -d [:space:])
						elif [ "${first}" = "D" ]
						then
							domain=$(echo "${line}" | cut -d'	' -f2 | tr -d [:space:])
						elif [ "${first}" = "()" ]
						then
							ssource=$(echo "${line}" | cut -d'(' -f2 | cut -d ')' -f1)
						fi
				done < "${main_dir}/${sample}/${sample}.tax"
				break
			elif [[ "${counter}" -gt 12 ]]; then
				echo "${main_dir}/${sample}/${sample}.tax does not exist, can not continue"
				exit 1
			else
				counter=$(( counter + 1 ))
				total=$(( counter * 5 ))
				echo "Waiting 5 seconds (total=${total}s) for ${main_dir}/${sample}/${sample}.tax to appear"
				sleep 5
			fi
		done
		# Check to see if prokka finished successfully

		# Set default busco database as bacteria in event that we dont have a database match for sample lineage
		buscoDB="bacteria_odb9"
		# Iterate through taxon levels (species to domain) and test if a match occurs to entry in database. If so, compare against it
		busco_found=0
		for tax in $species $genus $family $order $class $phylum $domain
		do
			if [ -d "${share}/DBs/BUSCO/${tax,}_odb9" ]
			then
				buscoDB="${tax}_odb9"
				busco_found=1
				break
			fi
		done
		# Report an unknown sample to the maintenance file to look into
		if [[ "${busco_found}" -eq 0 ]]; then
			lobal_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%ss")
			echo "BUSCO: ${domain} ${phylum} ${class} ${order} ${family} ${genus} ${species} - Found as ${project}/${sample} on ${global_time}" >> "${shareScript}/maintenance_To_Do.txt"
		fi
		# Show which database entry will be used for comparison
		echo "buscoDB:${buscoDB}"

		# Make script to run BUSCO
		echo -e "#!/bin/bash -l\n" > "${main_dir}/BUSCO_${sample}_${start_time}.sh"
		echo -e "#$ -o BUSCO_${sample}.out" >> "${main_dir}/BUSCO_${sample}_${start_time}.sh"
		echo -e "#$ -e BUSCO_${sample}.err" >> "${main_dir}/BUSCO_${sample}_${start_time}.sh"
		echo -e "#$ -N BUSCO_${sample}"   >> "${main_dir}/BUSCO_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/BUSCO_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/BUSCO_${sample}_${start_time}.sh"
		echo -e "module load busco/3.0.1\n" >> "${main_dir}/BUSCO_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/do_busco.sh\" \"${sample}\" ${buscoDB} \"${project}\"" >> "${main_dir}/BUSCO_${sample}_${start_time}.sh"

		mlst_db="false"
		while IFS= read -r line;
		do
			db_genus=$(echo "${line}" | cut -d "," -f1)
			db_species=$(echo "${line}" | cut -d',' -f2)
		 	echo "db_genus; ${db_genus}	:	db_species; ${db_species}	:	genus; ${genus}	:	species; ${species}"
		 	if [ "${db_species}" == "${species}" ] && [[ "${db_genus}" == "${genus}" ]]; then
		 		mlst_db="true"
				break
			fi
	 	done < "${share}/DBs/mlst_dbs.csv"

		if [[ "${mlst_db}" == "true" ]]; then
			# Need to check if database exists and the proper one if multiple
			# Create script to call SRST2 MLST
			if [[ "${genus}_${species}" == "Escherchia_coli" ]] || [[ "${genus}_${species}" == "Acinetobacter_baumannii" ]]; then
				temp_genspecies="${genus}_${species}#1"
				if [[ "${genus}" == "Escherchia" ]]; then
					alt_db_name="Achtman"
				else
					alt_db_name="Oxford"
				fi
				echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -o srst2MLST_${sample}.out" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -e srst2MLST_${sample}.err" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -N srst2MLST_${sample}"   >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -pe smp ${procs}" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "module load srst2/0.1.7\n" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "module unload Python/2.7.11" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "module load Python/2.7.15\n" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_srst2_mlst.sh\" \"${sample}\" \"${project}\" \"${genus}\" \"${species}#1\"" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "mv \"${main_dir}/${sample}/MLST/srst2/${sample}_srst2_${temp_genspecies}.mlst\" \"${main_dir}/${sample}/MLST/${sample}_srst2_${genus}_${species}_${alt_db_name}.mlst\"" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"

				temp_genspecies="${genus}_${species}#2"
				alt_db_name="Pasteur"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "#$ -o srst22MLST_${sample}.out" >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "#$ -e srst22MLST_${sample}.err" >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "#$ -N srst22MLST_${sample}"   >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "#$ -pe smp ${procs}" >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "module load srst2/0.1.7\n" >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "module unload Python/2.7.11" >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "module load Python/2.7.15\n" >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_srst2_mlst.sh\" \"${sample}\" \"${project}\" \"${genus}\" \"${species}#2\"" >> "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
				echo -e "mv \"${main_dir}/${sample}/MLST/srst2/${sample}_srst2_${temp_genspecies}.mlst\" \"${main_dir}/${sample}/MLST/${sample}_srst2_${genus}_${species}_${alt_db_name}.mlst\"" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
			else
				echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -o srst2MLST_${sample}.out" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -e srst2MLST_${sample}.err" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -N srst2MLST_${sample}"   >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -pe smp ${procs}" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "module load srst2/0.1.7\n" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "module unload Python/2.7.11" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "module load Python/2.7.15\n" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_srst2_mlst.sh\" \"${sample}\" \"${project}\" \"${genus}\" \"${species}\"" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
				echo -e "mv \"${main_dir}/${sample}/MLST/srst2/${sample}_srst2_${genus}_${species}.mlst\" \"${main_dir}/${sample}/MLST/${sample}_srst2_${genus}_${species}.mlst\"" >> "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
			fi
		else
			echo "${genus} ${species} does not have a pubMLST DB scheme"
		fi

		# Make script to run serotyper for E. coli only
		if [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
			echo -e "#!/bin/bash -l\n" > "${main_dir}/serty_${sample}_${start_time}.sh"
			echo -e "#$ -o serty_${sample}.out" >> "${main_dir}/serty_${sample}_${start_time}.sh"
			echo -e "#$ -e serty_${sample}.err" >> "${main_dir}/serty_${sample}_${start_time}.sh"
			echo -e "#$ -N serty_${sample}"   >> "${main_dir}/serty_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/serty_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/serty_${sample}_${start_time}.sh"
			# echo "module load serotyper\n" >> "${main_dir}/serty_${sample}_${start_time}.sh"
			echo -e "\"${shareScript}/run_serotyper.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/serty_${sample}_${start_time}.sh"
		fi
	done
}

# Loop 6 - Submit scripts for tools requiring species identification
submit_relies_on_species_confirmation() {
	check_time=$(date)
	echo "Loop 6 started at ${check_time}" >> "${processed}/${project}/${project}.log"
	for sample in "${sample_names[@]}";
	do
		echo "Loop 6 Send: ${sample}"
		prokka_success="false"
		counter=0
		total=0
		while :
		do
			if [[ -s "${main_dir}/${sample}/prokka/${sample}_PROKKA.gbf" ]] && [[ -s "${main_dir}/${sample}/prokka/${sample}_PROKKA.gff" ]]; then
				prokka_success="true"
				break
			elif [[ "${counter}" -gt 240 ]]; then
				break
			else
				counter=$(( counter + 1 ))
				total=$(( counter * 5 ))
				echo "Waiting 5 seconds (total=${total}s) for PROKKA output file(s) to appear"
				sleep 5
			fi
		done
		# Cleanup L5 scripts
		echo 'Cleaning L5-species determining scripts'
		mv "${main_dir}/taxID_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/taxID_${sample}_${start_time}.sh"
		mv "${shareScript}/taxID_${sample}"* "${main_dir}/${sample}/scripts/"
		echo -e 'Trying to call L6-species_dependant qsubs' #\n${main_dir}/BUSCO_${sample}_${start_time}.sh\n${main_dir}/srst2MLST_${sample}_${start_time}.sh\n(O)${main_dir}/serty_${sample}_${start_time}.sh'
		if [[ ${prokka_success} == "true" ]]; then
			qsub "${main_dir}/BUSCO_${sample}_${start_time}.sh"
		else
			echo "20 minutes has elapsed and ${main_dir}/${sample}/prokka.[gbf/gff] has not appeared"
		fi
		if [[ -f "${main_dir}/srst2MLST_${sample}_${start_time}.sh" ]]; then
			qsub "${main_dir}/srst2MLST_${sample}_${start_time}.sh"
		fi
		if [[ -f "${main_dir}/srst22MLST_${sample}_${start_time}.sh" ]]; then
			qsub "${main_dir}/srst22MLST_${sample}_${start_time}.sh"
		fi
		if [[ -f "${main_dir}/serty_${sample}_${start_time}.sh" ]]; then
			#qsub "${main_dir}/serty_${sample}_${start_time}.sh"
			:
		fi
	done
	check_time=$(date)
	echo "Loop 6 ended at ${check_time}" >> "${processed}/${project}/${project}.log"
}

# Loop 7 - Write scripts for validating run
make_validation_request() {
	for sample in "${sample_names[@]}";
	do
		echo "Loop 7 Make: ${sample}"

		# Make script to validate sample
		echo -e "#!/bin/bash -l\n" > "${main_dir}/validate_${sample}_${start_time}.sh"
		echo -e "#$ -o validate_${sample}.out" >> "${main_dir}/validate_${sample}_${start_time}.sh"
		echo -e "#$ -e validate_${sample}.err" >> "${main_dir}/validate_${sample}_${start_time}.sh"
		echo -e "#$ -N validate_${sample}"   >> "${main_dir}/validate_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/validate_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/validate_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/validate_piperun.sh\" \"${sample}\" \"${project}\" > \"${main_dir}/${sample}/${sample}_pipeline_stats.txt\"" >> "${main_dir}/validate_${sample}_${start_time}.sh"
	done
}

# Loop 7 - Submit scripts for validation
submit_validation_request() {
	check_time=$(date)
	echo "Loop 7 started at ${check_time}" >> "${processed}/${project}/${project}.log"
	for sample in "${sample_names[@]}";
	do
		echo "Loop 7 Send: ${sample}"
		validation_success="false"
		counter=0
		total=0
		if [[ "${genus}_${species}" == "Acinetobacter_baumannii" ]] || [[ "${genus}_${species}" == "Escherichia_coli" ]]; then
			if [[ "${species}" == "coli" ]]; then
				alt_db_name="Achtman"
			elif [[ "${species}" == "Oxford" ]]; then
				alt_db_name="Oxford"
			fi
			while :
			do
				if [[ -f "${main_dir}/${sample}/BUSCO/short_summary_${sample}.txt" ]] && [[ -f "${main_dir}/${sample}/MLST/${sample}_srst2_${genus}_${species}_Pasteur" ]] && [[ -f "${main_dir}/${sample}/MLST/${sample}_srst2_${genus}_${species}_${alt_db_name}" ]]; then
					Loop6_success="true"
					break
				elif [[ "${counter}" -gt 120 ]]; then
					break
				else
					counter=$(( counter + 1 ))
					total=$(( counter * 5 ))
					echo "Waiting 5 seconds (total=${total}s) for species dependant tools to complete (currently on ${sample})"
					sleep 5
				fi
			done
		elif [[ "${mlst_db}" == "true" ]]; then
			while :
			do
				if [[ -f "${main_dir}/${sample}/BUSCO/short_summary_${sample}.txt" ]] && [[ -f "${main_dir}/${sample}/MLST/${sample}_srst2_${genus}_${species}" ]]; then
					Loop6_success="true"
					break
				elif [[ "${counter}" -gt 120 ]]; then
					break
				else
					counter=$(( counter + 1 ))
					total=$(( counter * 5 ))
					echo "Waiting 5 seconds (total=${total}s) for species dependant tools to complete (currently on ${sample})"
					sleep 5
				fi
			done
		else
			while :
			do
				if [[ -f "${main_dir}/${sample}/BUSCO/short_summary_${sample}.txt" ]]; then
					Loop6_success="true"
					break
				elif [[ "${counter}" -gt 120 ]]; then
					break
				else
					counter=$(( counter + 1 ))
					total=$(( counter * 5 ))
					echo "Waiting 5 seconds (total=${total}s) for species dependant tools to complete (currently on ${sample})"
					sleep 5
				fi
			done
		fi
		if [[ "${Loop6_success}" == "true" ]]; then
			echo 'Cleaning L6-species_dependant scripts'
			mv "${main_dir}/BUSCO_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/BUSCO_${sample}_${start_time}.sh"
			mv "${shareScript}/BUSCO_${sample}"* "${main_dir}/${sample}/scripts/"
			mv "${main_dir}/srst2MLST_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/srst2MLST_${sample}_${start_time}.sh"
			mv "${shareScript}/srst2MLST_${sample}"* "${main_dir}/${sample}/scripts/"
			if [[ -f "${main_dir}/srst22MLST_${sample}_${start_time}.sh" ]]; then
				mv "${main_dir}/srst22MLST_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/srst22MLST_${sample}_${start_time}.sh"
				mv "${shareScript}/srst22MLST_${sample}"* "${main_dir}/${sample}/scripts/"
			fi
			if [[ -f "${main_dir}/serty_${sample}_${start_time}.sh" ]]; then
				mv "${main_dir}/serty_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/serty_${sample}_${start_time}.sh"
				mv "${shareScript}/serty_${sample}"* "${main_dir}/${sample}/scripts/"
			fi
			echo 'Trying to call L7-validation qsubs'
			qsub "${main_dir}/validate_${sample}_${start_time}.sh"
		else
			echo "Validation of ${sample} never finished"
		fi
	done
	check_time=$(date)
	echo "Loop 7 ended at ${check_time}" >> "${processed}/${project}/${project}.log"
}

# Check that all piperun stats are complete before running summary
summarize() {
	ready="false"
	for sample in "${sample_names[@]}";
	do
		counter=0
		total=0
		validation_success="false"
		pipeline_file="${main_dir}/${sample}/${sample}_pipeline_stats.txt"
		last_line=$(tail -n1 "${pipeline_file}" | cut -d' ' -f1)
		while :
		do
			if [[ "${last_line}" == "----------" ]]; then
				validation_success="true"
				break
			elif [[ "${counter}" -gt 240 ]]; then
				break
			else
				counter=$(( counter + 1 ))
				total=$(( counter * 5 ))
				echo "Waiting 5 seconds (total=${total}s) for all validations to complete (currently on ${sample})"
				sleep 5
				last_line=$(tail -n1 "${pipeline_file}" | cut -d' ' -f1)
				echo ":${last_line}:"
			fi
		done
		echo "Checking that all samples are finished validating"
		if [[ "${validation_success}" == "true" ]]; then
			echo -e 'Cleaning L7-validation scripts'
			mv "${main_dir}/validate_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/validate_${sample}_${start_time}.sh"
			mv "${shareScript}/validate_${sample}"* "${main_dir}/${sample}/scripts/"
			ready="true"
		else
			echo "${sample} did not finish validation in 20 minutes"
		fi
	done

	if [[ "${ready}" = "true" ]]; then
		echo -e "#!/bin/bash -l\n" > "${main_dir}/SUM_${sample}_${start_time}.sh"
		echo -e "#$ -o SUM_${sample}.out" >> "${main_dir}/SUM_${sample}_${start_time}.sh"
		echo -e "#$ -e SUM_${sample}.err" >> "${main_dir}/SUM_${sample}_${start_time}.sh"
		echo -e "#$ -N SUM_${sample}"   >> "${main_dir}/SUM_${sample}_${start_time}.sh"
		echo -e "#$ -cwd"  >> "${main_dir}/SUM_${sample}_${start_time}.sh"
		echo -e "#$ -q short.q\n"  >> "${main_dir}/SUM_${sample}_${start_time}.sh"
		echo -e "\"${shareScript}/run_sum.sh\" \"${project}\"" >> "${main_dir}/SUM_${sample}_${start_time}.sh"
		echo -e "echo $(date) > \"${main_dir}/${sample}/logs/${sample}_summary_complete.txt\"" >> "${main_dir}/SUM_${sample}_${start_time}.sh"

		qsub "${main_dir}/SUM_${sample}_${start_time}.sh"
	else
		echo "Run is not ready to be summarized"
	fi

	### Figure out how to pipe output
	counter=0
	total=0
	summary_success="false"
	while :
	do
		if [[ -f "${main_dir}/${sample}/logs/${sample}_summary_complete.txt" ]]; then
			summary_success="true"
			break
		elif [[ "${counter}" -gt 60 ]]; then
			break
		else
			counter=$(( counter + 1 ))
			total=$(( counter * 5 ))
			echo "Waiting 5 seconds (total=${total}s) for summary to complete"
			sleep 5
		fi
	done

	if [[ "${summary_success}" == "true" ]]; then
		runsum=$(${shareScript}/view_sum.sh ${project})
		mv "${main_dir}/SUM_${sample}_${start_time}.sh" "${main_dir}/${sample}/scripts/SUM_${sample}_${start_time}.sh"
		mv "${shareScript}/SUM_${sample}"* "${main_dir}/${sample}/scripts/"
		echo "${runsum}"
	else
		echo "Summary did not complete"
	fi
}

report_completion() {
	# Add print time the run completed in the text that will be emailed
	global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%ss")
	echo "Finished with run ${project} at ${global_end_time}"
	runsum+=("
	${project} finished at ${global_end_time}")

	#Send email to submitter and Nick with run status
	requestor=$(whoami)
	if [ "${requestor}" != "nvx4" ]; then
		echo "Sending summary email to ${requestor}@cdc.gov & nvx4@cdc.gov"
		printf "%s\\n" "${runsum[@]}" | mail -s "Run Status for ${project}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
		printf "%s\\n" "${runsum[@]}" | mail -s "Run Status for ${project}_on_${run_start_time}_run.log" "${requestor}@cdc.gov"
	else
		echo "Sending summary email to nvx4@cdc.gov"
		printf "%s\\n" "${runsum[@]}" | mail -s "Run Status for ${project}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
	fi
}

if [[ "${do_download}" == "true" ]]; then
	if [[ "${assembly_on}" == "true" ]]; then
		do_assembly_download
	else
		#do_reads_download
		:
	fi
fi

#make_list_from_folder
make_list_from_list
# Loop 1
#make_fastq_unzipper
#submit_fastq_unzipper
# Loop 2
#make_relies_on_unzipped_fastqs
#submit_relies_on_unzipped_fastqs
# Loop 3
#make_relies_on_trimmed_fastqs
#submit_relies_on_trimmed_fastqs

check_for_Assemblies
exit
# Loop 4

make_relies_on_trimmed_assemblies
submit_relies_on_trimmed_assemblies
# Loop 5
make_relies_on_species_files
submit_relies_on_species_files
# Loop 6
make_relies_on_species_confirmation
submit_relies_on_species_confirmation
# Loop 7
make_validation_request
submit_validation_request

#summarize
#report_completion

rm ${shareScript}/config.sh

exit
