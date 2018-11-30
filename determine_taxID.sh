#!/bin/sh -l

#$ -o getTax.out
#$ -e getTax.err
#$ -N getTax
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Creates a summary file for the run and prints out a one word status for each sample in the run
#
# Usage ./run_sum.sh run_id
#
# No modules required
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_sum.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample_id supplied to determine_taxID.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./determine_taxID.sh sample_ID run_ID"
	echo "Output is saved to ${processed}/run_ID/sample_ID/taxonomy.csv"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty run_id supplied to determine_taxID.sh, exiting"
	exit 1
fi

sample=${1}
project=${2}

Domain="Not assigned"
Phylum="Not assigned"
Class="Not assigned"
Order="Not assigned"
Family="Not assigned"
Genus="Not assigned"
species="Not assigned"

if [[ -s "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_All).txt" ]]; then
	source="ANI"
	#echo "${source}"
	# Lookup Taxonomy
  #echo "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_All).txt"
	header=$(head -n 1 "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_All).txt")
	#echo "${header}"
	Genus=$(echo "${header}" | cut -d' ' -f1 | cut -d'-' -f2)
	species=$(echo "${header}" | cut -d' ' -f2 | cut -d'(' -f1)
	echo "${Genus}-${species}"
elif [[ -s "${processed}/${project}/${sample}/kraken/postAssembly/${sample}_kraken_summary_assembled_BP_data.txt" ]]; then
	source="Kraken"
	#echo "${source}"
	while IFS= read -r line;
	do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "S" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{print $4}')
		elif [ "${first}" = "G" ]
		then
			Genus=$(echo "${line}" | awk -F ' ' '{print $4}')
		fi
	#	elif [ "${first}" = "F" ]
	#	then
	#		family=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "O" ]
	#	then
	#		order=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "C" ]
	#	then
	#		class=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "P" ]
	#	then
	#		phylum=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "K" ]
	#	then
	#		domain=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "D" ]
	#	then
	#		domain=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	fi
	done < "${processed}/${project}/${sample}/kraken/postAssembly/${sample}_kraken_summary_assembled_BP_data.txt"
elif [[ -s "${processed}/${project}/${sample}/gottcha/${sample}_gottcha_species_summary.txt" ]]; then
	source="GOTTCHA"
	#echo "${source}"
	while IFS= read -r line;
	do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "S" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{print $5}')
		elif [ "${first}" = "G" ]
		then
			Genus=$(echo "${line}" | awk -F ' ' '{print $4}')
		fi
	#	elif [ "${first}" = "F" ]
	#	then
	#		family=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "O" ]
	#	then
	#		order=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "C" ]
	#	then
	#		class=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "P" ]
	#	then
	#		phylum=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "K" ]
	#	then
	#		domain=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	elif [ "${first}" = "D" ]
	#	then
	#		domain=$(echo "${line}" | awk -F ' ' '{print $4}')
	#	fi
	done < "${processed}/${project}/${sample}/gottcha/${sample}_gottcha_species_summary.txt"
elif [[ -s "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt" ]]; then
	source="16s"
	#echo "${source}"
		# Lookup Taxonomy
		line=$(tail -n 1 "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt")
		Genus=$(echo "${line}" | cut -d"	" -f3 | cut -d" " -f1)
		species=$(echo "${line}" | cut -d"	" -f3 | cut -d" " -f2)
else
	echo "Exiting, no reliable taxonomy files exist"
	exit 1
fi

if [[ ${Genus} == "Peptoclostridium" ]]; then
	Genus="Clostridium"
fi

while IFS= read -r line;
do
	DB_genus=$(echo ${line} | cut -d"," -f1)
	#echo ":${Genus}:${DB_genus}:"
	if [[ "${Genus,}" = "${DB_genus}" ]]; then
			tax_DB="${local_DBs}/taxes.csv"
			Domain=$(echo "${line}" | cut -d"," -f2)
			Phylum=$(echo "${line}" | cut -d"," -f3)
			Class=$(echo "${line}" | cut -d"," -f4)
			Order=$(echo "${line}" | cut -d"," -f5)
			Family=$(echo "${line}" | cut -d"," -f6 | tr -d '\r' )
			#echo ":${Family}:"
			break
	fi
done < "${local_DBs}/taxes.csv"
printf "(${source}) \nD:	${Domain}\nP:	${Phylum}\nC:	${Class}\nO:	${Order}\nF:	${Family}\nG:	${Genus}\ns:	${species}\n" > "${processed}/${project}/${sample}/${sample}.tax"
