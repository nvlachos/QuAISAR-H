#!/bin/sh -l

#$ -o run_c-sstar.out
#$ -e run_c-sstar.err
#$ -N run_c-sstar
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
# No modules needed

#
# Finds anti-microbial resistance genes in the resFinder and ARG-ANNOT databases and exports a file containing list of all genes found
#
# Usage ./run_c-sstar.sh sample_name run_type(g/u for gapped/ungapped) similarity(l/m/h/u/p/o for low(80),medium(95),high(98),ultra-high(99),perfect(100),other(set in config.sh)) miseq_run_id plasmid(optional)
#
# No Modules needed
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_c-sstar.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_c-sstar.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_c-sstar.sh   sample_name   run-type[g/u](for gapped/ungapped)   similarity[l/m/h/u/p/o](for low/medium/high/ultra-high/perfect as 80/95/98/99/100, other(st in config.sh) miseq_run_id path/to/DB plasmid(optional)"
	echo "Output is saved to ${processed}/sample_name/c-sstar"
	exit
elif [ -z "$2" ]; then
	echo "Empty run type supplied to run_c-sstar.sh, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty similarity supplied to run_c-sstar.sh, exiting"
	exit 1
elif [ -z "$4" ]; then
	echo "Empty project id supplied to run_c-sstar.sh, exiting"
	exit 1
elif [ -z "$5" ] || [ ! -f "${5}" ]; then
	echo "Empty alternate ID supplied to run_c-sstar.sh, exiting"
	exit 1
fi

# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${4}/${1}"


#Set similarity threshold (set in config.sh) to values given in arguments
if [ "${3}" == "h" ]; then
	sim=${csstar_high}
elif [ "${3}" == "l" ]; then
	sim=${csstar_low}
elif [ "${3}" == "u" ]; then
	sim=${csstar_ultrahigh}
elif [ "${3}" == "m" ]; then
	sim=${csstar_medium}
elif [ "${3}" == "p" ]; then
	sim=${csstar_perfect}
elif [ "${3}" == "o" ]; then
	sim=${csstar_other}
else
	echo "Unknown similarity threshold set (use 'l,m,h,u,or p' for 80,95,98,99,or 100% respectively). Defaulting to 98%"
	sim=${csstar_high}
fi
if [[ "${6}" == "--plasmid" ]] || [[ "${6}" == "-p" ]]; then
	if [[ -s "${OUTDATADIR}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed.fasta" ]]; then
		source_assembly="${OUTDATADIR}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed.fasta"
		OUTDATADIR="${OUTDATADIR}/c-sstar_plasmid"
	else
		"No anti-microbial genes were found using c-SSTAR because there were No Plasmids Found" > "${OUTDATADIR}/c-sstar_plasmid/${1}_plasmid_scaffolds_trimmed.fasta"
		exit
	fi
else
	source_assembly="${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta"
	OUTDATADIR="${OUTDATADIR}/c-sstar"
fi

#alt_database="${5##*/}"
alt_database_path=$(basename -- "${5}")
alt_database=$(echo ${alt_database_path##*/} | cut -d'.' -f1)
alt_database=${alt_database//_srst2/}
echo ${alt_database}

# Creates the output c-sstar folder if it does not exist yet
#echo "${OUTDATADIR}"
if [ ! -d "$OUTDATADIR" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR"
	mkdir -p "$OUTDATADIR"
fi

# Set and call proper version of script based upon if gaps are allowed or not
if [ "${2}" == "u" ]; then
	suffix="ungapped"
	if [ ! -d "$OUTDATADIR/${alt_database}" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/${alt_database}"
		mkdir -p "$OUTDATADIR/${alt_database}_${suffix}"
	fi
	owd=$(pwd)
	cd "${OUTDATADIR}/${alt_database}_${suffix}"
	echo "Running c-SSTAR on ResGANNOT DB"
	python "${shareScript}/c-SSTAR_ungapped.py" -g "${source_assembly}" -s "${sim}" -d "${5}" > "${OUTDATADIR}/${alt_database}_${suffix}/${1}.${alt_database}.${suffix}_${sim}.sstar"
elif [ "${2}" == "g" ]; then
	suffix="gapped"
	if [ ! -d "$OUTDATADIR/${alt_database}_${suffix}" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/${alt_database}_${suffix}"
		mkdir -p "$OUTDATADIR/${alt_database}_${suffix}"
	fi
	owd=$(pwd)
	cd "${OUTDATADIR}/${alt_database}_${suffix}"
	echo "Running c-SSTAR on ResGANNOT DB"
	python "${shareScript}/c-SSTAR_gapped.py" -g "${source_assembly}" -s "${sim}" -d "${5}" > "${OUTDATADIR}/${alt_database}_${suffix}/${1}.${alt_database}.${suffix}_${sim}.sstar"
else
	echo "Unknown run type set (only use 'g' or 'u' for gapped/ungapped analysis"
	exit 1
fi



###################################### FIND WAY TO CATCH FAILURE !!!!!!!!!! ###############################

# Goes through ResGANNOT outfile and adds labels as well as resistance conferred to the beginning of the line
# Takes .sstar file in and outputs as .sstar_grouped
while IFS= read -r line
do
	line=${line,,}
	#echo ${line}
	label1=$(echo "${line}" | cut -d '	' -f3 | tr '[:upper:]' '[:lower:]')
	label2=$(echo "${line}" | cut -d '	' -f4 | tr '[:upper:]' '[:lower:]')
	info1=""
	#info2=""
	#echo "R1;${label1}-${label2}"
	if [[ "${label1}" = *"tr" ]] && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 was truncated"
		label1="${label1:0:${#label1} - 2}"
		info1="${info1}trunc-"
	fi
	if ( [[ "${label1}" = *"*"* ]] || [[ "${label1}" = *"*" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is likely novel"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}novel-"
	fi
	if ( [[ "${label1}" = *"?"* ]] || [[ "${label1}" = *"?" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is uncertain due to incomplete alignment"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}alinc-"
	fi
	if ( [[ "${label1}" = *"$"* ]] || [[ "${label1}" = *"$" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is uncertain due to incomplete alignment"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}edge-"
	fi
	label1=$(echo "${label1,,}" | tr -d '*?$')
	label2=$(echo "${label2,,}" | tr -d '*?$')
	source=$(echo "${line}" | cut -d '	' -f1 | tr -d '[:space:]')
	resistance=$(echo "${line}" | cut -d '	' -f2 | tr -d '[:space:]')
	contig=$(echo "${line}" | cut -d '	' -f5 | tr -d '[:space:]')
	percent=$(echo "${line}" | cut -d '	' -f6 | cut -d'%' -f1 | tr -d '[:space:]')
	len1=$(echo "${line}" | cut -d '	' -f7 | tr -d '[:space:]')
	len2=$(echo "${line}" | cut -d '	' -f8 | tr -d '[:space:]')
	SNPs=$(echo "${line}" | cut -d '	' -f10 | tr -d '[:space:]')
	plen=$(echo "${line}" | cut -d '	' -f9 | tr -d '[:space:]')
	#if [[ ${len1} -ge ${len2} ]]; then
	#	plen=100
	#else
	#	plen=$( echo "${len1} ${len2}" | awk '{ printf "%2.f", ($1*100)/$2 }' )
	#fi
	if [[ -z "${info1}" ]]; then
		info1="normal"
	else
		info1=${info1::-1}
	fi
	#printf "%-10s %-50s %-15s %-25s %-25s %-40s %-4s %-5d %-5d %-5d\\n" "${source}" "${resistance}" "${label1}" "${info1}" "${label2}" "${contig}" "${percent}" "${len1}" "${len2}" "${plen}" "${SNPs}"
	echo "${source}	${resistance}	${label1}	${info1}	${label2}	${contig}	${percent}	${len1}	${len2}	${plen}" "${SNPs}"
done < "${OUTDATADIR}/${alt_database}_${suffix}/${1}.${alt_database}.${suffix}_${sim}.sstar" > "${OUTDATADIR}/${alt_database}_${suffix}/${1}.${alt_database}.${suffix}_${sim}.sstar_grouped"
#echo "${OUTDATADIR}/${alt_database}_${suffix}/${1}.${alt_database}.${suffix}_${sim}.sstar_grouped"
#echo "${OUTDATADIR}/${1}.${alt_database}.${suffix}_${sim}_sstar_summary.txt"
sort -k7,7nr -k10,10nr -k8,8n "${OUTDATADIR}/${alt_database}_${suffix}/${1}.${alt_database}.${suffix}_${sim}.sstar_grouped" > "${OUTDATADIR}/${1}.${alt_database}.${suffix}_${sim}_sstar_summary.txt"
#sort -k7,7nr -k10,10nr -k8,8n "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/151217_M02103_0050_000000000-AKVY1/2013732586/c-sstar/ResGANNOT_20180608_gapped/2013732586.ResGANNOT_20180608.gapped_98.sstar_grouped" > "/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles/151217_M02103_0050_000000000-AKVY1/2013732586/c-sstar/2013732586.ResGANNOT_20180608.gapped_98_sstar_summary.txt"

# Catches an empty or missing file
if [ ! -s "${OUTDATADIR}/${1}.${alt_database}.${suffix}_${sim}_sstar_summary.txt" ]; then
	echo "No anti-microbial genes were found using c-SSTAR with both resFinder and ARG-ANNOT DBs" > "${OUTDATADIR}/${1}.${alt_database}.${suffix}_${sim}_sstar_summary.txt"
fi

#Returns to original directory
cd "${owd}"

#Script exited gracefully (unless something else inside failed)
exit 0
