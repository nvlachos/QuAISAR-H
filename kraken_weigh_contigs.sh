#!/bin/sh -l

#$ -o kraken_weigh_contigs.out
#$ -e kraken_weigh_contigs.err
#$ -N kraken_weigh_contigs
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
# No MODs needed

#
# Grabs the best species match based on %/read hits from the kraken tool run
#
# Usage ./kraken_weigh_contigs.sh sample_name run_id kraken_version[kraken|kraken2]
#
# No modules needed to run
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to kraken_weigh_contigs.sh, exiting"
	exit 1
elif [ -z "$1" ]; then
	echo "Empty sample name supplied to kraken_weigh_contigs.sh, exiting"
	exit 1
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./kraken_weigh_contigs.sh  sample_name  [pre/post] [paired/assembled] run_id"
	echo "Output is saved to ${processed}/miseq_run_id_id/sample_name/kraken/(pre/post)assembly/sample_name_kraken_summary_(paired/assembled)"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id supplied to kraken_weigh_contigs.sh, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty kraken version supplied to kraken_weigh_contigs.sh, exiting"
	exit 1
fi

#Sets output folder to the correct path relative to assembly completion
OUTDATADIR="${processed}/${2}/${1}/${3}/postAssembly"
echo "-${OUTDATADIR}-"


#Checks to see if the list file used for calculations exists and exits if it does not
if [[ ! -s "${OUTDATADIR}/${1}_assembled_BP.${3}" ]]; then
	echo "${OUTDATADIR}/${1}_assembled_BP.${3} does not exist"
	exit 1
fi

if [[ "${3}" != "kraken" ]] && [[ "${3}" != "kraken2" ]]; then
	echo "Invalid kraken version supplied, exiting"
	exit 4
fi

contig_sizes=[]
total_size=0
unclassified=0

#Parses the kraken output list line by line
while IFS= read -r line
do
		classified==$(echo "${line}" | cut -d'	' -f1)
		if [[ "${classified}" == "C" ]]; then
			contig_size=$(echo "${line}" | cut -d'	' -f4)
			contig_sizes+=()${contig_size})
			total_size=$(( total_size + contig_size ))
		else
			echo "Contig not classified"
			unclassified=$(( unclassified + 1 ))
		fi
done < "${OUTDATADIR}/${1}_assembled_BP.${3}"

contig_count=${#contig_sizes[@]}
echo "Contig count = ${contig_count}"
echo "Total Size = ${total_size}"
echo "unclassified = ${unclassified}"

if [[ ! -s "${OUTDATADIR}/${1}_assembled_weighted.mpa" ]]; then
	echo "${OUTDATADIR}/${1}_assembled_weighted.mpa does not exist, cant do mpa adjustment"
	exit 1
fi





#Script exited gracefully (unless something else inside failed)
exit 0
