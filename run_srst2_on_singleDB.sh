#!/bin/sh -l

#$ -o srst2.out
#$ -e srst2.err
#$ -N srst2
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/prep_srst2.sh"

#
# Usage ./run_srst2.sh   sample_name   run_ID
#
# script uses srst2 to find AR genes from ResGANNOT DBs.
#
#  Modules needed (loaded via prep_srst2.sh script - Python/2.7.15, bowtie2/2.2.9, samtools/0.1.18
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_srtst2.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./run_srst2.sh  sample_name run_ID"
	echo "Output location is ${processed}/run_ID/sample_name/srst2"
	exit 0
fi

# Create output directory
mkdir "${processed}/${2}/${1}/srst2"

# Check if the proper files exist to perform srst2
# Find R1 zipped file, or zip the raw trimmed reads, and copy to output directory
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" ]; then
		#echo "1"
		cp "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" ]; then
		#echo "2"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz"
	elif [[ ! -d "${processed}/${2}/${1}/trimmed" ]]; then
		#echo "5"
		if [[ -f "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" ]] && [[ ! -f "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq" ]]; then
			gunzip -c "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" > "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq"
		fi
		if [[ -f "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" ]] && [[ ! -f "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq" ]]; then
			gunzip -c "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" > "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq"
		fi
		bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq" in2="${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq" out="${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" out2="${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
		trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" "${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R1_001.unpaired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
		cat "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}.paired.fq"
		cat "${processed}/${2}/${1}/trimmed/${1}_R1_001.unpaired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.unpaired.fq" > "${processed}/${2}/${1}/trimmed/${1}.single.fq"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	fi
fi
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" ]; then
		#echo "3"
		cp "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" ]; then
		#echo "4"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz"
	fi
fi

# Prints the command that will be submitted to use srst2 to find AR genes
echo "--input_pe ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz --output ${processed}/${2}/${1}/srst2/${1}_ResGANNOT --gene_db ${resGANNOT_srst2}"

# Calls srst2 with the options for AR discovery
#python2 "${shareScript}/srst2/scripts/srst2.py" --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_ResGANNOT" --gene_db "${resGANNOT_srst2}"
srst2 --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_ResGANNOT" --gene_db "${resGANNOT_srst2}"

# Cleans up leftover files
#rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
#rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
rm -r "${processed}/${2}/${1}/srst2/"*".bam"
rm -r "${processed}/${2}/${1}/srst2/"*".pileup"

# Removes the extra ResGANNOT__ from all files created
find ${processed}/${2}/${1}/srst2 -type f -name "*ResGANNOT__*" | while read FILE ; do
  dirname=$(dirname $FILE)
	filename=$(basename $FILE)
	filename="${filename/_ResGANNOT__/__}"
	#echo "Found-${FILE}"
	#echo "${filename}"
    mv "${FILE}" "${dirname}/${filename}"
done

# Close out modules that were loaded specifically for srst2
. "${mod_changers}/close_srst2.sh"
