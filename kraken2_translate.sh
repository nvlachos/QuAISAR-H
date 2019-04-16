#!/bin/sh -l

#$ -o kraken2_translate.out
#$ -e best_hit_from_kraken.err
#$ -N best_hit_from_kraken
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
# Usage ./best_hit_from_kraken.sh sample_name pre/post(relative to assembly) source_type(paired/assembled) run_id [full]
#
# No modules needed to run
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to best_hit_from_kraken.sh, exiting"
	exit 1
elif [ -z "$1" ]; then
	echo "Empty sample name supplied to best_hit_from_kraken.sh, exiting"
	exit 1
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./best_hit_from_kraken.sh  sample_name  [pre/post] [paired/assembled] run_id"
	echo "Output is saved to ${processed}/miseq_run_id_id/sample_name/kraken/(pre/post)assembly/sample_name_kraken_summary_(paired/assembled)"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id supplied to best_hit_from_kraken.sh, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty assembly relativity supplied to best_hit_from_kraken.sh, exiting"
	exit 1
fi

#Sets output folder to the correct path relative to assembly completion
OUTDATADIR="${processed}/${2}/${1}/kraken2/${3}Assembly"
echo "-${OUTDATADIR}-"
if [[ "${3}" == "post" ]]; then
	post_fix="assembled"
elif [[ "${3}" == "pre" ]]; then
	post_fix="paired"
else
	echo "Unknown relativity ${3}, exiting"
fi

#Checks to see if the list file used for calculations exists and exits if it does not
if [[ ! -s "${OUTDATADIR}/${1}_${post_fix}.kraken2" ]]; then
	echo "${OUTDATADIR}/${1}_${post_fix}.list does not exist"
	exit 1
fi

> ${OUTDATADIR}/${1}_${post_fix}.labels

who_am_i=$(whoami)

# Create associative array to hold taxid and tree styructure so not to look up as many times
declare -a tax_trees

#Parses the kraken output list line by line
while IFS= read -r line
do
		contig_info=$(echo "${line}" | cut -d'	' -f2)
		taxid=$(echo "${line}" | cut -d'	' -f3)
		#echo "${contig_info} ::::: ${taxid}"



		if [[ -n ${tax_trees[$taxonomy]} ]]; then
			echo "Looking up" 
			taxonomy=$(python ${shareScript}/entrez_get_taxon_from_number.py ${taxid} ${who_am_i})
			tax_tree=$(echo "${taxonomy}" | cut -d'|' -f3 | cut -d'	' -f2)
		else
			print "Already in array"
			tax_tree=${tax_trees[${taxonomy}]}
		fi
		echo "${contig_info} ::: ${tax_tree}"

done < "${OUTDATADIR}/${1}_${post_fix}.kraken2"

#Script exited gracefully (unless something else inside failed)
exit 0
