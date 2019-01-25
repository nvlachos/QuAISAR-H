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
# Creates a single file that attempts to pull the best taxonomic information from the isolate. Currently, it operates in a linera fashion, e.g. 1.ANI, 2.kraken, 3.Gottcha, 4.16s
# The taxon is chosen based on the highest ranked classifier first
#
# Usage ./determine_texID.sh sample_name project_ID
#
# No modules required
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied todetermine_taxID.sh, exiting"
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

# Set default values after setting variables
sample=${1}
project=${2}

Domain="Not_assigned"
Phylum="Not_assigned"
Class="Not_assigned"
Order="Not_assigned"
Family="Not_assigned"
Genus="Not_assigned"
species="Not_assigned"
source="Not_assigned"
confidence_index="0"
source_file="Not_assigned"


Check_source() {
	start_at="${1}"
	if [[ "${start_at}" -le 1 ]]; then
		for f in ${processed}/${project}/${sample}/ANI/*; do
			if [[ "${f}" = *"best_ANI_hits_ordered"* ]]; then
				header=$(head -n1 ${f})
				if [[ ${header} != "No matching ANI database found for"* ]]; then
		    	do_ANI
		    	return
				fi
			fi
		done
	fi
	if [[ "${start_at}" -le 2 ]]; then
		if [[ -s "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt" ]]; then
			best_line=$(head -n1 "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt")
			largest_line=$(tail -n1 "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt")
			IFS='	' read -r -a best_array <<< "$best_line"
			IFS='	' read -r -a largest_array <<< "$largest_line"
			best_arr_size="${#best_array[@]}"
			largest_arr_size="${#largest_array[@]}"
			best_species=$(echo ${best_line} | cut -d'	' -f3)
			largest_species=$(echo ${largest_line} | cut -d'	' -f3)
			#echo "largest:${largest_species}:"
			#echo "best:${best_species}:"
			if [[ "${largest_arr_size}" -ge 3 ]]; then
				do_16s "largest"
				return
			elif [[ "${best_arr_size}" -ge 3 ]]; then
				do_16s "best"
				return
			fi
		fi
	fi
	if [[ "${start_at}" -le 3 ]];then
		if [[ -s "${processed}/${project}/${sample}/gottcha/${sample}_gottcha_species_summary.txt" ]]; then
			do_GOTTCHA
			return
		fi
	fi
	if [[ "${start_at}" -le 4 ]]; then
		if [[ -s "${processed}/${project}/${sample}/kraken/postAssembly/${sample}_kraken_summary_assembled_BP_data.txt" ]]; then
			do_Kraken
		return
		fi
	fi
	echo "No ACCEPTABLE taxonomy file found to use"
}

do_ANI() {
	source="ANI"
	#echo "${source}"
	if [[ -f "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_All).txt" ]]; then
		source_file="${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered(${sample}_vs_All).txt"
	else
		source_file=$(ls -t "${processed}/${project}/${sample}/ANI/best_ANI_hits_ordered"* | head -n 1)
	fi
	header=$(head -n 1 "${source_file}")
	#echo "${header}"
	Genus=$(echo "${header}" | cut -d' ' -f1 | cut -d'-' -f2)
	species=$(echo "${header}" | cut -d' ' -f2 | cut -d'(' -f1)
	confidence_index=$(echo "${header}" | cut -d' ' -f1 | cut -d'%' -f1)
	#echo "${Genus}-${species}"
}

do_16s() {
	if [[ "${1}" = "largest" ]]; then
		source_file="${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt"
		line=$(tail -n 1 "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt")
		source="16s_largest"
		confidence_index=$(head -n1 "${processed}/${project}/${sample}/16s/${sample}.nt.RemoteBLASTN.sorted" | cut -d'	' -f3)
		confidence_index="${confidence_index}"
	elif [[ "${1}" = "best" ]]; then
		line=$(head -n 1 "${processed}/${project}/${sample}/16s/${sample}_16s_blast_id.txt")
		source="16s_best"
		confidence_index=$(head -n1 "${processed}/${project}/${sample}/16s/${sample}.nt.RemoteBLASTN" | cut -d'	' -f3)
		confidence_index="${confidence_index}"
	else
		break
	fi
	Genus=$(echo "${line}" | cut -d"	" -f3 | cut -d" " -f1)
	species=$(echo "${line}" | cut -d"	" -f3 | cut -d" " -f2)
}

do_GOTTCHA() {
	source="GOTTCHA"
	source_file="${processed}/${project}/${sample}/gottcha/${sample}_gottcha_species_summary.txt"
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
		confidence_index=$(tail -n1 "${source_file}" | cut -d' ' -f2)
		confidence_index="${confidence_index}"
	done < "${processed}/${project}/${sample}/gottcha/${sample}_gottcha_species_summary.txt"
}

do_Kraken() {
	source="Kraken"
	source_file="${processed}/${project}/${sample}/kraken/postAssembly/${sample}_kraken_summary_assembled_BP_data.txt"
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
	done < "${processed}/${project}/${sample}/kraken/postAssembly/${sample}_kraken_summary_assembled_BP_data.txt"
	confidence_index=$(tail -n1 "${source_file}" | cut -d' ' -f2)
	confidence_index="${confidence_index}"
}

Check_source 0

if [[ ! -z ${Genus} ]]; then
	Genus=$(echo ${Genus} | tr -d [:space:])
fi
if [[ ! -z ${species} ]]; then
	species=$(echo ${species} | tr -d [:space:])
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

printf "(${source})-${confidence_index}%%-${source_file}\nD:	${Domain}\nP:	${Phylum}\nC:	${Class}\nO:	${Order}\nF:	${Family}\nG:	${Genus}\ns:	${species}\n" > "${processed}/${project}/${sample}/${sample}.tax"
