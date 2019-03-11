#!/bin/sh -l

#$ -o get_run_plasFlow.out
#$ -e get_run_plasFlow.err
#$ -N get_run_plasFlow
#$ -cwd
#$ -q all.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#. "${mod_changers}/pipeline_mods"
#. ./module_changers/pipeline_mods
#. ./module_changers/list_modules.sh

#
# Will attempt to find any plasmids in sample
#
# Usage ./run_plasFlow.sh sample_name run_id
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_plasmidFinder.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_plasmidFinder.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_plasmFlow.sh  sample_name run_id"
	echo "Output by default is ${processed}/miseq_run_id/sample_name/plasmFlow"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty run_id supplied to run_plasFlow.sh, exiting"
	exit 1
fi

module load PlasFlow/1.1
module load Python/3.5.4
module load bowtie2/2.2.9
module load samtools/1.4.1
module load bam2fastq/1.1.0
module load Unicycler/0.4.4;
module load gcc/5.5;
module load SPAdes/3.11.1;
module load racon/1.2.0;


. ./module_changers/list_modules.sh

# Create output directory
if [[ ! -d "${processed}/${2}/${1}/plasFlow" ]]; then
	mkdir "${processed}/${2}/${1}/plasFlow"
fi

if [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" ]; then
	echo "1"
	:
elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" ]; then
	echo "2"
	gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz"
else
	echo "3"
	if [[ -f "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" ]] && [[ ! -f "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq" ]]; then
			gunzip -c "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" > "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq"
	fi
	if [[ -f "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" ]] && [[ ! -f "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq" ]]; then
		gunzip -c "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" > "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq"
	fi
	echo "Running BBDUK and trimmomatic"
	bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq" in2="${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq" out="${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" out2="${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
	trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" "${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R1_001.unpaired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
fi



# Trim contigs a little to 2000 and larger and put through plasflow.
# The remaining analysis steps have not yet been completed
if [[ -s "${processed}/${2}/${1}/Assembly/${1}_scaffolds_trimmed.fasta" ]]; then
	python2 "${shareScript}/removeShortContigs.py" "${processed}/${2}/${1}/Assembly/scaffolds.fasta" "2000"
	${shareScript}/fasta_headers.py "${processed}/${2}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${processed}/${2}/${1}/plasFlow/${1}_scaffolds_trimmed_2000.fasta"
	rm -r "${processed}/${2}/${1}/Assembly/scaffolds.fasta.TRIMMED.fasta"
	PlasFlow.py --input "${processed}/${2}/${1}/plasFlow/${1}_scaffolds_trimmed_2000.fasta" --output "${processed}/${2}/${1}/plasFlow/${1}_plasFlow_results.tsv" --threshold 0.7
	mkdir ${processed}/${2}/${1}/plasFlow/bowtie2-index/
	bowtie2-build -f "${processed}/${2}/${1}/plasFlow/${1}_plasFlow_results.tsv_chromosomes.fasta" "${processed}/${2}/${1}/plasFlow/bowtie2-index/bowtie2_${1}_chr"
	${processed}/${2}/${1}/plasFlow/filtered_reads_70/
	bowtie2 -x "${processed}/${2}/${1}/plasFlow/bowtie2-index/bowtie2_${1}_chr" -1 "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" 2 "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" -S "${processed}/${2}/${1}/plasFlow/filtered_reads_70/${1}.sam" -p 12 --local
	samtools view -bS "${processed}/${2}/${1}/plasFlow/filtered_reads_70/${1}.sam" > "${processed}/${2}/${1}/plasFlow/filtered_reads_70/${1}.bam"
	bam2fastq --no-aligned -o "${processed}/${2}/${1}/plasFlow/filtered_reads_70/${1}_R#_bacterial.fastq" "${processed}/${2}/${1}/plasFlow/filtered_reads_70/${1}.bam"
	unicycler -1 "${processed}/${2}/${1}/plasFlow/filtered_reads_70/${1}_R_1_bacterial.fastq" -2 "${processed}/${2}/${1}/plasFlow/filtered_reads_70/${1}_R_2_bacterial.fastq" -o "${processed}/${2}/${1}/plasFlow/Unicycler_assemblies/${1}_uni_assembly"
fi

module unload PlasFlow/1.1
module unload Python/3.5.4
module unload bowtie2/2.2.9
module unload samtools/1.4.1
module unload bam2fastq/1.1.0
module unload Unicycler/0.4.4;
module unload gcc/5.5;
module unload SPAdes/3.11.1;
module unload racon/1.2.0;
