#!/bin/bash -l

#$ -o primary_processing.out
#$ -e primary_processing.err
#$ -N ppX
#$ -pe smp 14
#$ -cwd
#$ -q all.q

# Copy config file into working directory to allow changes to made to output directory if necessary
shareScript=$(pwd)
echo "${shareScript}"
if [[ -f "${shareScript}/config_template.sh" ]]; then
	if [[ -f "${shareScript}/config.sh" ]]; then
		rm -r "${shareScript}/config.sh"
	fi
	echo "Trying to copy config_template.sh"
	cp "${shareScript}/config_template.sh" "${shareScript}/config.sh"
fi

#Import the config file with shortcuts and settings
. ${shareScript}/config.sh

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
			is_full_run="true"
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
			if [[ -z ${INDATADIR} ]]; then
				INDATADIR="${processed}/${PROJECT}/"
			fi
			if [[ ${BASEDIR} = "${requestor}_${global_time}" ]]; then
				BASEDIR="${processed}"
			fi
			list_path="${processed}/${PROJECT}/${PROJECT}_list.txt"
			run_name="project_${PROJECT}"
			do_download="true"
			is_full_run="true"
			shift 2
			;;
		#Tells the script to run analyses from already downloaded fastq files
		-n | --no_download)
			do_download=false
			shift
			;;
		#Tells the script that only the files found in the attached list need to be run
		-l | --list)
			list_path="$2"
			quick_list=$(echo "${2}" | cut -d'.' -f1)
			#INDATADIR="${processed}/${PROJECT}" # NOT USED yet in list mode
			if [[ -z ${BASEDIR} ]]; then
				BASEDIR="${processed}"
			fi
			do_download="false"
			run_name="list_${quick_list}"
			is_full_run="false"
			shift 2
			;;
		#Tells the script that only the single isolate needs to be run
		-s | --single)
			echo "${2}/${3}" > "./tempList_${global_time}.txt"
			list_path="./tempList_${global_time}.txt"
			INDATADIR="${processed}/${3}"
			if [[ -z ${BASEDIR} ]]; then
				BASEDIR="${processed}"
			fi
			trn=$(echo "${2}" | sed 's/\//-/')
			run_name="single_${trn}"
			do_download="false"
			is_full_run="false"
			shift 3
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
echo -e "Source folder: ${INDATADIR}\\nOutput folder: ${BASEDIR}\\nDownload fastqs(.gzs): ${do_download}\\nList based analysis:  ${list_path}"

# Checks that a full FASTQ source path is given
if [ "${INDATADIR:0:1}" != "/" ] && [ "${INDATADIR:0:1}" != "$" ]; then
	echo "${INDATADIR}"
	echo 'ERROR: The full path was not specified.' >&2
	exit 1
fi

# Actual Pipeline of tools and analyses
process_samples()	{

	#Time tracker to gauge time used by each step
	totaltime=0

	# Set arguments to filename(sample name) project (miseq run id) and outdatadir(${processed}/project/filename)
	filename="${1}"
	project="${2}"
	OUTDATADIR="${3}/${2}"

	# Create an empty time_summary file that tracks clock time of tools used
	touch "${OUTDATADIR}/${filename}/${filename}_time_summary.txt"
	time_summary=${OUTDATADIR}/${filename}/${filename}_time_summary.txt

	echo "Time summary for ${project}/${filename}:" >> "${time_summary}"
	echo "${project}/${filename} started at ${global_time}"

	echo "Starting processing of ${project}/${filename}"
	#Checks if FASTQ folder exists for current sample
	if [[ -d "$OUTDATADIR/$filename/FASTQs" ]]; then
		# Checks if FASTQ folder contains any files then continue
		if [[ "$(ls -A "${OUTDATADIR}/${filename}/FASTQs")" ]]; then
			# Checks to see if those files in the folder are unzipped fastqs
			count_unzip=`ls -1 ${OUTDATADIR}/${filename}/FASTQs/*.fastq 2>/dev/null | wc -l`
			count_zip=`ls -1 ${OUTDATADIR}/${filename}/FASTQs/*.fastq.gz 2>/dev/null | wc -l`
			if [[ ${count_unzip} != 0 ]]; then
			#if [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}"*".fastq" ]]; then
				echo "----- FASTQ(s) exist, continuing analysis -----"
				if [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq" ]] && [[ ! -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq.gz" ]]; then
					gzip < "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq.gz"
				fi
				if [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq" ]] && [[ ! -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" ]]; then
					gzip < "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz"
				fi
			# Checks if they are zipped fastqs (checks for R1 first)
			elif [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq.gz" ]]; then
				#echo "R1 zipped exists - unzipping"
				gunzip -c "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq.gz" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq"
				# Checks for paired R2 file
				if [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" ]]; then
					#echo "R2 zipped exists - unzipping"
					gunzip -c "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq"
				else
					echo "No matching R2 to unzip :("
				fi
			# Checks to see if there is an abandoned R2 zipped fastq
			elif [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" ]]; then
				#echo "R2 zipped  exists - unzipping"
				gunzip -c "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq"
				echo "No matching R1 to unzip :("
			fi
		# If the folder is empty then return from function
		else
			echo "FASTQs folder empty - No fastqs available for ${filename} (and download was not requested). Either unzip fastqs to $OUTDATADIR/FASTQs or run the -d flag to trigger unzipping of gzs"
			return 1
		fi
	# If the fastq folder does not exist then return out of function
	else
		echo "FASTQs not downloaded and FASTQs folder does not exist for ${filename}. No fastqs available (and download was not requested). Unzip fastqs to ${OUTDATADIR}/FASTQs"
		return 1
	fi

	# Get start time for qc check
	start=$SECONDS
	### Count the number of Q20, Q30, bases and reads within a pair of FASTQ files
	echo "----- Counting read quality -----"
	# Checks for and creates the specified output folder for the QC counts
	if [ ! -d "$OUTDATADIR/$filename/preQCcounts" ]; then
		echo "Creating $OUTDATADIR/$filename/preQCcounts"
		mkdir -p "$OUTDATADIR/$filename/preQCcounts"
	fi
	# Run qc count check on raw reads
	python "${shareScript}/Fastq_Quality_Printer.py" "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq" "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq" > "${OUTDATADIR}/$filename/preQCcounts/${filename}_counts.txt"

	# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeQCcount=$((end - start))
	echo "QC count - ${timeQCcount} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeQCcount))

	###  Trimming and Quality Control  ###
	echo "----- Running BBDUK on reads -----"
	# Gets start time for bbduk
	start=$SECONDS
	# Creates folder for BBDUK output
	if [ ! -d "$OUTDATADIR/$filename/removedAdapters" ]; then
		echo "Creating $OUTDATADIR/$filename/removedAdapters"
		mkdir -p "$OUTDATADIR/$filename/removedAdapters"
	# It complains if a folder already exists, so the current one is removed (shouldnt happen anymore as each analysis will move old runs to new folder)
	else
		echo "Removing old $OUTDATADIR/$filename/removedAdapters"
		rm -r "$OUTDATADIR/$filename/removedAdapters"
		echo "Recreating $OUTDATADIR/$filename/removedAdapters"
		mkdir -p "$OUTDATADIR/$filename/removedAdapters"
	fi
	# Run bbduk
	bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq" in2="${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq" out="${OUTDATADIR}/${filename}/removedAdapters/${filename}-noPhiX-R1.fsq" out2="${OUTDATADIR}/${filename}/removedAdapters/${filename}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
	# Get end time of bbduk and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeAdapt=$((end - start))
	echo "Removing Adapters - ${timeAdapt} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeAdapt))

	### Quality and Adapter Trimming using trimmomatic ###
	echo "----- Running Trimmomatic on reads -----"
	# Get start time of trimmomatic
	start=$SECONDS
	# Creates folder for trimmomatic output if it does not exist
	if [ ! -d "$OUTDATADIR/$filename/trimmed" ]; then
		mkdir -p "$OUTDATADIR/$filename/trimmed"
	fi
	# Run trimmomatic
	trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${OUTDATADIR}/${filename}/removedAdapters/${filename}-noPhiX-R1.fsq" "${OUTDATADIR}/${filename}/removedAdapters/${filename}-noPhiX-R2.fsq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.paired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.unpaired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.paired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
	# Get end time of trimmomatic and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeTrim=$((end - start))
	echo "Trimming - ${timeTrim} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeTrim))


	# Check differences after QC and trimming (also for gottcha proper read count for assessing unclassified reads)
	# Get start time for qc check on trimmed reads
	start=$SECONDS
	### Count the number of Q20, Q30, bases and reads within the trimmed pair of FASTQ files
	echo "----- Counting read quality of trimmed files-----"
	# Checks for and creates the specified output folder for the QC counts
	if [ ! -d "$OUTDATADIR/$filename/preQCcounts" ]; then
		echo "Creating $OUTDATADIR/$filename/preQCcounts"
		mkdir -p "$OUTDATADIR/$filename/preQCcounts"
	fi
	# Run qc count check on filtered reads
	python "${shareScript}/Fastq_Quality_Printer.py" "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.paired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.paired.fq" > "${OUTDATADIR}/${filename}/preQCcounts/${filename}_trimmed_counts.txt"

		# Merge both unpaired fq files into one for GOTTCHA
	cat "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.paired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.paired.fq" > "${OUTDATADIR}/${filename}/trimmed/${filename}.paired.fq"
	cat "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.unpaired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.unpaired.fq" > "${OUTDATADIR}/${filename}/trimmed/${filename}.single.fq"
	#gzip < "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.paired.fq" > "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.paired.fq.gz"
	#gzip < "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.paired.fq" > "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.paired.fq.gz"

		# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeQCcount_trimmed=$((end - start))
	echo "QC count trimmed - ${timeQCcount_trimmed} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeQCcount))



	######  Run Kraken on cleaned reads  ######
	echo "----- Running Kraken on cleaned reads -----"
	# Get start time of kraken on reads
	start=$SECONDS
	# Run kraken
	"${shareScript}/run_kraken.sh" "${filename}" pre paired "${project}"
	# Get end time of kraken on reads and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeKrak=$((end - start))
	echo "Kraken - ${timeKrak} seconds" >> "${time_summary}"
 	#${shareScript}/run_kraken.sh ${filename} pre R1
 	#${shareScript}/run_kraken.sh ${filename} pre R2
 	#${shareScript}/run_kraken.sh ${filename} pre single
	totaltime=$((totaltime + timeKrak))

	##### Run gottcha(v1) on cleaned reads #####
	echo "----- Running gottcha on cleaned reads -----"
	# Get start time of gottcha
	start=$SECONDS
	# run gootcha
	"${shareScript}/run_gottcha.sh" "${filename}" "${project}"
	# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeGott=$((end - start))
	echo "Gottcha - ${timeGott} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeGott))

	# Check reads using SRST2
	echo "----- Running SRST2 -----"
	start=$SECONDS
	"${shareScript}/run_srst2_on_singleDB.sh" "${filename}" "${project}"
	"${shareScript}/run_srst2_on_singleDB_alternateDB.sh" "${filename}" "${project}" "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"
	end=$SECONDS
	timesrst2=$((end - start))
	echo "SRST2 - ${timesrst2} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timesrst2))

	######  Assembling Using SPAdes  ######
	echo "----- Assembling Using SPAdes -----"
	# Get start time of SPAdes
	start=$SECONDS
	# script tries 3 times for a completed assembly
	for i in 1 2 3
	do
		#Deletes any core dumps from SPAdes failing

		# If assembly exists already and this is the first attempt (then the previous run will be used) [should not happen anymore as old runs are now renamed]
		if [ -s "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" ]; then
			echo "Previous assembly already exists, using it (delete/rename the assembly folder at ${OUTDATADIR}/ if you'd like to try to reassemble"
		# Run normal mode if no assembly file was found
		else
			"${shareScript}/run_SPAdes.sh" "${filename}" normal "${project}"
		fi
		# Removes any core dump files (Occured often during testing and tweaking of memory parameter
		if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
			echo "Found core dump files in assembly (assumed to be from SPAdes, but could be anything before that as well) and attempting to delete"
			find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
		fi
	done
	# Returns if all 3 assembly attempts fail
	if [[ -f "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" ]] && [[ -s "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" ]]; then
		echo "Assembly completed and created a non-empty scaffolds file"
	else
		echo "Assembly FAILED 3 times, continuing to next sample..." >&2
		return 1
	fi
	for i in 1 2 3
	do
		# If plasmid Assembly exists already and this is the first attempt (then the previous run will be used) [should not happen anymore as old runs are now renamed]
		if [ -f "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta" ]; then
			echo "Previous plasmid assembly already exists, using it (delete/rename the assembly folder at ${OUTDATADIR}/ if you'd like to try to reassemble"
		else
			"${shareScript}/run_SPAdes.sh" "${filename}" plasmid "${project}"
			# Removes any core dump files created from SPAdes (occurred fairly regularly during testing/tweaking)
			if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
				echo "Found core dump files in plasmid Assembly (assumed to be from SPAdes, but could be anything before that as well) and attempting to delete"
				find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;

			else
				echo "No core dump files during plasmid assembly attempt number ${i}, plasmidSPAdes finished successfully and found nothing, creating dummy scaffolds file"
				>> "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta"
			fi
		fi
	done
	# Returns if all 3 assembly attempts fail
	if [ ! -f "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta" ]; then
		echo "plasmid Assembly FAILED 3 times, continuing to next step..." >&2
	fi
	# Get end time of SPAdes and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeSPAdes=$((end - start))
	echo "SPAdes - ${timeSPAdes} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeSPAdes))

	### Removing Short Contigs  ###
	echo "----- Removing Short Contigs -----"
	python "${shareScript}/removeShortContigs.py" "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" 500
	mv "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta"

	### Removing Short Contigs  ###
	echo "----- Removing Short plasmid Contigs -----"
	python "${shareScript}/removeShortContigs.py" "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta" 500
	mv "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta" "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed.fasta"

	# Checks to see that the trimming and renaming worked properly, returns if unsuccessful
	if [ ! -s "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta" ]; then
		echo "Trimmed contigs file does not exist continuing to next sample">&2
		return 1
	fi

	### ReKraken on Assembly ###
	echo "----- Running Kraken on Assembly -----"
	# Get start time of kraken on assembly
	start=$SECONDS
	# Run kraken on assembly
	"${shareScript}/run_kraken.sh" "${filename}" post assembled "${project}"
	# Get end time of kraken on assembly and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeKrakAss=$((end - start))
	echo "Kraken on Assembly - ${timeKrakAss} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeKrakAss))

	# Get ID fom 16s
	echo "----- Identifying via 16s blast -----"
	start=$SECONDS
	"${shareScript}/16s_blast.sh" "${filename}" "${project}"
	end=$SECONDS
	time16s=$((end - start))
	echo "16S - ${time16s} seconds" >> "${time_summary}"
	totaltime=$((totaltime + time16s))

	# Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
	"${shareScript}/determine_taxID.sh" "${filename}" "${project}"

	# Capture the anticipated taxonomy of the sample using kraken on assembly output
	echo "----- Extracting Taxonomy from Taxon Summary -----"
	# Checks to see if the kraken on assembly completed successfully
	if [ -s "${OUTDATADIR}/${filename}/${filename}.tax" ]; then
		# Read each line of the kraken summary file and pull out each level  taxonomic unit and store for use later in busco and ANI
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
			elif [ "${first}" = "F" ]
			then
				family=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "O" ]
			then
				order=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "C" ]
			then
				class=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "P" ]
			then
				phylum=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "K" ]
			then
				kingdom=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "D" ]
			then
				domain=$(echo "${line}" | awk -F ' ' '{print $2}')
			fi
		done < "${OUTDATADIR}/${filename}/${filename}.tax"
	# Print out taxonomy for confirmation/fun
	echo "Taxonomy - ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species}"
	# If no kraken summary file was found
	else
		echo "No Taxonomy output available to make best call from, skipped"
	fi

	### Check quality of Assembly ###
	echo "----- Running quality checks on Assembly -----"
	# Get start time of QC assembly check
	start=$SECONDS
	# Run qc assembly check
	"${shareScript}/run_Assembly_Quality_Check.sh" "${filename}" "${project}"
	# Get end time of qc quality check and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeQCcheck=$((end - start))
	echo "QC Check of Assembly - ${timeQCcheck} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeQCcheck))

	### Prokka on assembly ###
	echo "----- Running Prokka on Assembly -----"
	# Get start time for prokka
	start=$SECONDS
	# Run prokka
	"${shareScript}/run_prokka.sh" "${filename}" "${project}"
	# Get end time of prokka and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeProk=$((end - start))
	echo "Identifying Genes (Prokka) - ${timeProk} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeProk))

	# Rename contigs to something helpful (Had to wait until after prokka runs due to the strict naming requirements
	mv "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed_original.fasta"
	python3 "${shareScript}/fasta_headers.py" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed_original.fasta" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta"
	if [[ -s "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed.fasta" ]]; then
		mv "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed.fasta" "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed_original.fasta"
		python3 "${shareScript}/fasta_headers.py" "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed_original.fasta" "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed.fasta"
	fi

	### Average Nucleotide Identity ###
	echo "----- Running ANI for Species confirmation -----"
	# ANI uses assembly and sample would have exited already if assembly did not complete, so no need to check
	# Get start time of ANI
	start=$SECONDS
	# run ANI
	# Temp fix for strange genera until we do vs ALL all the time.
	if [[ "${genus}" = "Peptoclostridium" ]] || [[ "${genus}" = "Clostridioides" ]]; then
		genus="Clostridium"
	fi
	"${shareScript}/run_ANI.sh" "${filename}" "${genus}" "${species}" "${project}"
	#"${shareScript}/run_ANI.sh" "${filename}" "All" "All" "${project}"
	# Get end time of ANI and calculate run time and append to time summary (and sum to total time used
	end=$SECONDS
	timeANI=$((end - start))
	echo "autoANI - ${timeANI} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeANI))

	# Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
	"${shareScript}/determine_taxID.sh" "${filename}" "${project}"
	"${OUTDATADIR}/${filename}/${filename}.tax"

	### BUSCO on prokka output ###
	echo "----- Running BUSCO on Assembly -----"
	# Check to see if prokka finished successfully
	if [ -s "${OUTDATADIR}/${filename}/prokka/${filename}_PROKKA.gbf" ] || [ -s "${OUTDATADIR}/${filename}/prokka/${filename}_PROKKA.gff" ]; then
		# Get start time of busco
		start=$SECONDS
		# Set default busco database as bacteria in event that we dont have a database match for sample lineage
		buscoDB="bacteria_odb9"
		# Iterate through taxon levels (species to domain) and test if a match occurs to entry in database. If so, compare against it
		busco_found=0
		for tax in $species $genus $family $order $class $phylum $kingdom $domain
		do
			if [ -d "${local_DBs}/BUSCO/${tax,}_odb9" ]
			then
				buscoDB="${tax,}_odb9"
				busco_found=1
				break
			fi
		done
		# Report an unknown sample to the maintenance file to look into
		if [[ "${busco_found}" -eq 0 ]]; then
			global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
			echo "BUSCO: ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species} - Found as ${project}/${filename} on ${global_time}" >> "${shareScript}/maintenance_To_Do.txt"
		fi
		# Show which database entry will be used for comparison
		echo "buscoDB:${buscoDB}"
		# Run busco
		"${shareScript}/do_busco.sh" "${filename}" "${buscoDB}" "${project}"
		# Get end time of busco and calculate run time and append to time summary (and sum to total time used
		end=$SECONDS
		timeBUSCO=$((end - start))
		echo "BUSCO - ${timeBUSCO} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeBUSCO))
	# Prokka did not complete successfully and busco cant run (since it relies on prokka output)
	else
		echo "Prokka output not found, not able to process BUSCO"
	fi


	### c-SSTAR for finding AR Genes ###
	echo "----- Running c-SSTAR for AR Gene identification -----"
	# c-SSTAR uses assembly and sample would have exited already if assembly did not complete, so no need to check
	# Get start time of ccstar
	start=$SECONDS

	# Run csstar in default mode from config.sh
	"${shareScript}/run_c-sstar_on_single.sh" "${filename}" "${csstar_gapping}" "${csstar_identity}" "${project}"
	"${shareScript}/run_c-sstar_on_single_alternate_DB.sh" "${filename}" "${csstar_gapping}" "${csstar_identity}" "${project}" "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"
	# Should the parameters be different when checking on plasmids specifically
	"${shareScript}/run_c-sstar_on_single.sh" "${filename}" "${csstar_gapping}" "${csstar_plasmid_identity}" "${project}" "--plasmid"
	"${shareScript}/run_c-sstar_on_single_alternate_DB.sh" "${filename}" "${csstar_gapping}" "${csstar_plasmid_identity}" "${project}" "--plasmid" "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"

	# Get end time of csstar and calculate run time and append to time summary (and sum to total time used
	end=$SECONDS
	timestar=$((end - start))
	echo "c-SSTAR - ${timestar} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timestar))

	# Get MLST profile
	echo "----- Running MLST -----"
	start=$SECONDS
	"${shareScript}/run_MLST.sh" "${filename}" "${project}"
	if [[ "${genus}_${species}" = "Acinetobacter_baumannii" ]]; then
		"${shareScript}/run_MLST.sh" "${filename}" "${project}" "-f" "abaumannii"
	elif [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
		# Verify that ecoli_2 is default and change accordingly
		"${shareScript}/run_MLST.sh" "${filename}" "${project}" "-f" "ecoli_2"
	fi
	end=$SECONDS
	timeMLST=$((end - start))
	echo "MLST - ${timeMLST} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeMLST))

	# Try to find any plasmids
	echo "----- Identifying plasmids using plasmidFinder -----"
	start=$SECONDS
	"${shareScript}/run_plasmidFinder.sh" "${filename}" "${project}" "plasmid"
	"${shareScript}/run_plasmidFinder.sh" "${filename}" "${project}" "plasmid_on_plasmidAssembly"
	end=$SECONDS
	timeplasfin=$((end - start))
	echo "plasmidFinder - ${timeplasfin} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeplasfin))

	# Append total time to bottom of time summary
	echo "Total time: ${totaltime} seconds" >> "${time_summary}"

	# Designate end of this sample #
	global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
	echo "

				End of sample ${filename}
				completed at ${global_end_time}
	"

	# Extra dump cleanse in case anything else failed
	if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
		echo "Found core dump files at end of processing ${filename} and attempting to delete"
		find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
	fi

} # End process_samples

#Copies over and unzips all 'Determined' zipped FASTQ files from the sequencing instrument for the requested project_id
if [[ "${do_download}" = "true" ]]; then
	echo "Attempting to download fastq files..."
	if [[ "${indir_set}" = "true" ]]; then
		"${shareScript}/get_Reads_from_folder.sh" "${PROJECT}" "${INDATADIR}" "${postfix}"
	else
		bash "${shareScript}/get_Reads_from_Instrument.sh" "${PROJECT}"
	fi
else
	if [[ "${is_full_run}" = "true" ]]; then
		echo "Attempting to make list of files in MiSeqAnalysisFiles folder ${PROJECT}"
		# If a list does not exist, a temporary one with all files in folder is used to test which ones are folders that contain FASTQs sub-folders
		if [[ ! -f "${processed}/${PROJECT}/${PROJECT}_list.txt" ]]; then
			ls ${processed}/${PROJECT}/ > "${processed}/${PROJECT}/${PROJECT}_list_temp.txt"
			while IFS=' ' read -r line;
			do
				#echo "${line}"
				if [[ -d "${processed}/${PROJECT}/${line}/FASTQs" ]]; then
					echo "${PROJECT}/${line}" >> "${processed}/${PROJECT}/${PROJECT}_list.txt"
				else
					echo "Sample '${line}' does not contain a FASTQs folder"
				fi
			done < ${processed}/${PROJECT}/${PROJECT}_list_temp.txt
			rm "${processed}/${PROJECT}/${PROJECT}_list_temp.txt"
		else
			echo "List file already exists for ${PROJECT}"
		fi
	fi
fi

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
if [[ "${is_full_run}" = "true" ]]; then
	mkdir "${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}"
	log_file="${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}/${run_name}_on_${run_start_time}.log"
else
	mkdir "${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}"
	log_file="${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}/${run_name}_on_${run_start_time}.log"
fi

#Get the time the run started to use as the identifier
outarray=()
echo "Full Run started at ${run_start_time}   status file will be ${run_name}_on_${run_start_time}.log"
outarray+=("${run_name} started at ${run_start_time} and saved to ${run_name}_on_${run_start_time}_full_run.log")


#Each file in the list is put through the full pipeline
for projfile in "${file_list[@]}";
do
	#echo "${projfile}"
	file=$(echo "${projfile}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	proj=$(echo "${projfile}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	process_samples "${file}" "${proj}" "${BASEDIR}"
	"${shareScript}/validate_piperun.sh" "${file}" "${proj}" > "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"
	cat "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" >> "${log_file}"
	#echo "Done"
done

#Full run overview status file is created for quick view
# The full run log is searched and brief status overview is returned for each sample in the log as SUCCESSFUL,WARNING, or FAILED along with counts of each
if [[ "${is_full_run}" = "true" ]];then
	cp "${log_file}" "${processed}/${proj}/"
	runsumdate=$(date "+%m_%d_%Y_at_%Hh_%Mm")
	mv "${processed}/${PROJECT}/${run_name}_on_${run_start_time}.log" "${processed}/${PROJECT}/${PROJECT}_run_summary_at_${runsumdate}.sum"
	# summary file is consolidated and prepped to send to email recipient(s)
	runsum=$(${shareScript}/view_sum.sh ${PROJECT})
	outarray+="${runsum}"
# If a list or single sample was run through the pipeline instead of a full run. The log file is copied over and renamed
elif [[ "${is_full_run}" = "false" ]];then
	echo "In not full run move"
	runsumdate=$(date "+%m_%d_%Y_at_%Hh_%Mm")
	mv "${log_file}" "${Quaisar_H_log_directory}/${quick_list}_summary_at_${runsumdate}.sum"
	# summary file is consolidated and prepped to send to email recipient(s)
	runsum=$(${shareScript}/view_sum.sh ${quick_list}_summary_at_${runsumdate}.sum -l)
	outarray+="${runsum}"
fi

#Move log files and remove other temp files
# If job was submitted move out and err files to processing logs folder

if [[ ${host} = "cluster"* ]]; then
	qh_iteration=$(basename $0 | cut -d'.' -f1 | cut -d'_' -f3)
	if [[ -z "${qh_iteration}" ]]; then
		outbase="primary_processing"
	else
		outbase="primary_processing_${qh_iteration}"
	fi
	if [[ ! -d "${Quaisar_H_log_directory}" ]]; then
		mkdir "${Quaisar_H_log_directory}"
	fi
	if [ -f "${shareScript}/${outbase}.out" ]; then
		mv "${shareScript}/${outbase}.out" "${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}/"
		rename "${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}/${outbase}.out" "${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}/${run_name}_on_${run_start_time}_primary_processing.out"
	fi
	if [ -f "${shareScript}/${outbase}.err" ]; then
		mv "${shareScript}/${outbase}.err" "${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}/"
		rename "${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}/${outbase}.err" "${Quaisar_H_log_directory}/${run_name}_on_${run_start_time}/${run_name}_on_${run_start_time}_primary_processing.err"
	fi
fi


# Run the Seqlog creator on the proper file
if [ "${is_full_run}" = "true" ]; then
	"${shareScript}/make_Seqlog_from_log.sh" "${PROJECT}"
else
	if [ -s "${Quaisar_H_log_directory}/tempList.txt" ]; then
		"${shareScript}/make_Seqlog_from_list.sh" "${Quaisar_H_log_directory}/tempList.txt"
		rm "${Quaisar_H_log_directory}/tempList.txt"
	else
		"${shareScript}/make_Seqlog_from_list.sh" "${2}"
	fi
fi

# Clean all extraneous files and folders from samples that were run to save space
for projfile in "${file_list[@]}";
do
	#echo ${file}
	file=$(echo "${projfile}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	proj=$(echo "${projfile}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	"${shareScript}/sample_cleaner.sh" "${file}" "${proj}"
done

# Add print time the run completed in the text that will be emailed
global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Finished with run ${run_name} at ${global_end_time}"
outarray+=("
${run_name} finished at ${global_end_time}")
exit
#Send email to submitter and Nick with run status

if [ "${requestor}" != "nvx4" ]; then
	echo "Sending summary email to ${requestor}@cdc.gov & nvx4@cdc.gov"
	printf "%s\\n" "${outarray[@]}" | mail -s "Run Status for ${run_name}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
	printf "%s\\n" "${outarray[@]}" | mail -s "Run Status for ${run_name}_on_${run_start_time}_run.log" "${requestor}@cdc.gov"
else
	echo "Sending summary email to nvx4@cdc.gov"
	printf "%s\\n" "${outarray[@]}" | mail -s "Run Status for ${run_name}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
fi

# One final check for any dump files
if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
	echo "Found core dump files at end of processing ${run_name} and attempting to delete"
	find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
fi

if [[ -f "./config.sh" ]]; then
	#rm -r "./config.sh"
	:
fi
