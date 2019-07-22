#!/bin/sh -l

#$ -o srst2.out
#$ -e srst2.err
#$ -N srst2
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#Import the module file that loads all necessary mods
#. "${mod_changers}/prep_srst2.sh"

ml Python2/2.7.12 samtools/0.1.18 perl/5.16.1-MT srst2 bowtie2/2.2.4

#
# Usage ./run_srst2_on_singleDB_alternateDB.sh.sh   sample_name   MiSeq_Run_ID
#
# script uses srst2 to find AR genes from a custom srst2 formatted DB.
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./run_srst2_on_singleDB_alternateDB.sh.sh  sample_name MiSeq_Run_ID path_to_alt_DB"
	echo "Output location is ${processed}/run_ID/srst2"
	exit 0
fi

alt_DB_path=${3}
alt_DB=$(echo ${alt_DB_path##*/} | cut -d'.' -f1)
alt_DB=${alt_DB//_srst2/}

# Create output folder
if [[ ! -d "${processed}/${2}/${1}/srst2" ]]; then
	mkdir "${processed}/${2}/${1}/srst2"
fi

# Prep any reads that have not been trimmed
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" ]; then
		#echo "1"
		cp "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" ]; then
		#echo "2"
		$(gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz")
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz"
	fi
else
	echo "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz already exists"
fi
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" ]; then
		#echo "3"
		cp "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" ]; then
		#echo "4"
		$(gzip -c "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz")
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz"
	fi
else
	echo "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz already exists"
fi

# Display command that will be submitted
echo "--input_pe ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz --output ${processed}/${2}/${1}/srst2/${alt_DB} --gene_db ${alt_DB_path}"

# Submit command to run srst2 on alternate DB
python2 ${shareScript}/srst2-master/scripts/srst2.py --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_${alt_DB}" --gene_db "${alt_DB_path}"
#srst2 --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_${alt_DB}" --gene_db "${alt_DB_path}"

# Clean up extra files
rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
rm -r "${processed}/${2}/${1}/srst2/"*".bam"
rm -r "${processed}/${2}/${1}/srst2/"*".pileup"

# Rename files to reduce redundancy
find ${processed}/${2}/${1}/srst2 -type f -name "*_${alt_DB}__*" | while read FILE ; do
	#echo "Doing - $FILE"
  dirname=$(dirname $FILE)
	filename=$(basename $FILE)
	filename="${filename/_${alt_DB}__/__}"
	#echo "${filename}"
  mv "${FILE}" "${dirname}/${filename}"
done

#. "${mod_changers}/close_srst2.sh"
