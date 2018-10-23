#!/bin/sh -l

#$ -o srst2.out
#$ -e srst2.err
#$ -N srst2
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/scripts/config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"
#module clear

. "${mod_changers}/prep_srst2.sh"
. "${mod_changers}/list_modules.sh"


#
# Usage ./run_srst2.sh   sample_name   MiSeq_Run_ID
#
# script uses srst2 to find AR genes from resFinder and ARGANNOT DBs.
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_srtst2.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./run_srst2.sh  sample_name MiSeq_Run_ID Genus species"
	echo "Output location is ${processed}/${2}/srst2"
	exit 0
fi

if [[ ! -d "${processed}/${2}/${1}/srst2" ]]; then
	mkdir "${processed}/${2}/${1}/srst2"
fi

if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" ]; then
		#echo "1"
		cp "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" ]; then
		#echo "2"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz"
	elif [[ ! -d "${processed}/${2}/${1}/trimmed" ]]; then
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

if [ ! -d "${processed}/${2}/${1}/MLST/" ]; then
	mkdir "${processed}/${2}/${1}/MLST/"
	mkdir "${processed}/${2}/${1}/MLST/srst2"
elif [ ! -d "${processed}/${2}/${1}/MLST/srst2" ]; then
	mkdir "${processed}/${2}/${1}/MLST/srst2"
fi

cd "${processed}/${2}/${1}/MLST/srst2"

echo "do"
python2 "${shareScript}/srst2/scripts/getmlst.py" --species "${3} ${4}" > "${processed}/${2}/${1}/MLST/srst2/getmlst.out"
#getmlst.py --species "${3} ${4}" > "${processed}/${2}/${1}/MLST/srst2/getmlst.out"

echo "done"
if [[ "${3}" == "Acinetobacter" ]]; then
	echo "${processed}/${2}/${1}/MLST/srst2/${3}_${4}.fasta"
	if [[ "${4}" == "baumannii#1" ]]; then
		sed -i -e 's/Oxf_//g' "${processed}/${2}/${1}/MLST/srst2/${3}_${4}.fasta"
		sed -i -e 's/Oxf_//g' "${processed}/${2}/${1}/MLST/srst2/abaumannii.txt"
	elif [[ "${4}" == "baumannii#2" ]]; then
		sed -i -e 's/Pas_//g' "${processed}/${2}/${1}/MLST/srst2/${3}_${4}.fasta"
		sed -i -e 's/Pas_//g' "${processed}/${2}/${1}/MLST/srst2/abaumannii_2.txt"
	else
		echo "Unknown species in Acinetobacter MLST lookup"
	fi
fi

suggested_command=$(tail -n2 "${processed}/${2}/${1}/MLST/srst2/getmlst.out" | head -n1)
mlst_db=$(echo "${suggested_command}" | cut -d' ' -f11)
mlst_defs=$(echo "${suggested_command}" | cut -d' ' -f13)
mlst_delimiter=$(echo "${suggested_command}" | cut -d' ' -f15)
#echo "${mlst_db}"
#echo "${mlst_defs}"
#echo "${mlst_delimiter}"

if [[ "${mlst_delimiter}" != "'_'" ]]; then
	echo "Unknown delimiter - \"${mlst_delimiter}\""
	exit
else
	mlst_delimiter="_"
	#echo "Delimiter is OK (${mlst_delimiter})"
fi

echo "--input_pe ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz --output ${processed}/${2}/${1}/MLST/srst2 --mlst_db ${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}"

python "${shareScript}/srst2/scripts/srst2.py" --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/MLST/srst2/${1}" --mlst_db "${mlst_db}" --mlst_definitions "${mlst_defs}" --mlst_delimiter "${mlst_delimiter}"

mv "${processed}/${2}/${1}/MLST/srst2/${1}__mlst__${3}_${4}__results.txt" "${processed}/${2}/${1}/MLST/${1}_srst2_${3}_${4}.mlst"
rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
#rm -r "${processed}/${2}/${1}/MLST/srst2/*.tfa"
#rm -r "${processed}/${2}/${1}/MLST/${1}__${1}"*"scores"
if [[ -f "${processed}/${2}/${1}/MLST/srst2/${1}__${1}.${3}_${4}.pileup" ]]; then
	rm -r "${processed}/${2}/${1}/MLST/srst2/${1}__${1}.${3}_${4}.pileup"
fi
if [[ -f "${processed}/${2}/${1}/MLST/${1}__${1}.${3}_${4}.sorted.bam" ]]; then
	rm -r "${processed}/${2}/${1}/MLST/${1}__${1}.${3}_${4}.sorted.bam"
fi


. "${mod_changers}/close_srst2.sh"
