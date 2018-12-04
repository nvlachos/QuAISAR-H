#!/bin/sh -l

#$ -o run_csstar_proj_parser.out
#$ -e run_csstar_proj_parser.err
#$ -N run_csstar_proj_parser
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ -f config_template.sh ]]; then
	if [[ ! -f config.sh ]]; then
		cp config_template.sh config.sh
	fi
fi
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"
#. "${mod_changers}/list_modules.sh"


#
# Usage ./run_csstar_proj_parser.sh list_name gapped/ungapped (analysis ran) identity (80/95/98/99/100) output_name plasmid_identity(optional)
#
# Pulls out MLST, AR genes, and plasmid for the listed samples and consolidates them into one sheet (per category)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_csstar_proj_parser.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./run_csstar_proj_parser.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also) gapped/ungapped 80/95/98/99/100 output_prefix output_directory plasmid_identity_cutoff(optional, default = 40)"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${share}"
	exit 0
elif [[ ! -f ${1} ]]; then
	echo "list does not exit...exiting"
	exit 1
fi

output_directory=${5}
if [[ ! -d ${output_directory} ]]; then
	mkdir -p ${output_directory}
fi

# # Remove any pre-existing files from previous runs
if [[ -f ${output_directory}/${4}-mlst_summary.txt ]]; then
	rm ${output_directory}/${4}-mlst_summary.txt
fi
if [[ -f ${output_directory}/${4}-csstar_summary.txt ]]; then
	rm ${output_directory}/${4}-csstar_summary.txt
fi
if [[ -f ${output_directory}/${4}-plasmid_summary.txt ]]; then
	rm ${output_directory}/${4}-plasmid_summary.txt
fi
if [[ -f ${output_directory}/${4}_AR_plasmid_report.csv ]]; then
	rm ${output_directory}/${4}_AR_plasmid_report.csv
fi
if [[ -f ${output_directory}/${4}-csstar_summary_full.txt ]]; then
	rm ${output_directory}/${4}-csstar_summary_full.txt
fi
if [[ -f ${output_directory}/${4}-srst2.txt ]]; then
	rm ${output_directory}/${4}-srst2.txt
fi
if [[ -f ${output_directory}/${4}-srst2_rejects.txt ]]; then
	rm ${output_directory}/${4}-srst2_rejects.txt
fi

if [ "${3}" == 98 ]; then
	sim="h"
elif [ "${3}" == 80 ]; then
	sim="l"
elif [ "${3}" == 99 ]; then
	sim="u"
elif [ "${3}" == 95 ]; then
	sim="m"
elif [ "${3}" == 100 ]; then
	sim="p"
elif [ "${3}" == 40 ]; then
	sim="o"
fi
if [[ -z "${6}" ]]; then
	plaid=40
else
	plaid="${6}"
fi

rename="true"


declare -A groups
echo "Creating reference array"
counter=0
while IFS= read -r line;
do
	line=${line,,}
	gene=$(echo "${line}" | cut -d ':' -f1)
	first=${gene:0:1}
	if [ "$first" == "#" ]; then
		continue
	fi
	confers=$(echo "${line}" | cut -d ':' -f2)
	groups[${gene}]="${confers}"
	#echo "${counter}:${gene}:${confers}"
	counter=$(( counter + 1))
done < "${share}/DBs/star/group_defs.txt"

# Loop through and act on each sample name in the passed/provided list
 echo -e "\nUsing AR Database - ${resGANNOT_srst2_filename}\n"
 while IFS= read -r line; do
	 sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	 project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	 OUTDATADIR="${processed}/${project}/${sample_name}"
	 #rm -r "${processed}/${project}/${sample_name}/c-sstar/${resGANNOT_srst2_filename}"
	 #echo "Checking for ${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt"
	 if [[ -s ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq ]] && [[ ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq ]] || [[ -s ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz ]] && [[ ${processed}/${project}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
		 if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${resGANNOT_srst2_filename}_srst2__results.txt" ]] || [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${resGANNOT_srst2_filename}_srst2__results.txt" ]]; then
			 :
		 else
			 echo "It thinks it needs to put ${sample_name} through srst2"
			 #"${shareScript}/run_srst2_on_singleDB.sh" "${sample_name}" "${project}"
		 fi
	 fi

	 if [[ -s ${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt" ]];
		then
			:
			#echo "ResGANNOT file already exists for assembly of ${project}/${sample_name}"
		else
			#echo "${project}/${sample_name} - ${resGANNOT_srst_filename} - Not found"
			gapping=${2}
			gapping=${gapping:0:1}
			echo -e "Doing ResGANNOT as ${shareScript}/run_c-sstar_on_single.sh ${sample_name} ${gapping} ${sim} ${project}"
			"${shareScript}/run_c-sstar_on_single.sh" "${sample_name}" "${gapping}" "${sim}" "${project}"
		fi
	fi
	if [[ -s ${processed}/${project}/${sample_name}/plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta ]]; then
		#echo "Checking for - ${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${plaid}_sstar_summary.txt"
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${plaid}_sstar_summary.txt" ]];
		then
			:
			#echo "ResGANNOT file already exists for plasmidAssembly of ${project}/${sample_name}"
		else
			gapping=${2}
			gapping=${gapping:0:1}
			if [ "${plaid}" -eq 98 ]; then
				run_plaid="h"
			elif [ "${plaid}" -eq 80 ]; then
				run_plaid="l"
			elif [ "${plaid}" -eq 99 ]; then
				run_plaid="u"
			elif [ "${plaid}" -eq 95 ]; then
				run_plaid="m"
			elif [ "${plaid}" -eq 100 ]; then
				run_plaid="p"
			elif [ "${plaid}" -eq 40 ]; then
				run_plaid="o"
			fi
			echo -e "Doing ResGANNOT on plasmidAssembly as \n${shareScript}/run_c-sstar_on_single.sh ${sample_name} ${gapping} ${run_plaid} ${project} --plasmid"
			"${shareScript}/run_c-sstar_on_single.sh" "${sample_name}" "${gapping}" "${run_plaid}" "${project}" "--plasmid"
		fi
	else
		#echo "${project}/${sample_name} has no plasmid Assembly"
		:
	fi

	#ls ${processed}/${project}/${sample_name}/c-sstar_plasmid/

	sample_index=0
	oar_list=""
	# Looks at all the genes found for a sample
	if [[ -f "${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt" ]]; then
		ARDB_full="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
	else
		${shareScript}/run_c-sstar_on_single.sh "${sample_name}" "${gapping}" "${sim}" "${project}"
		ARDB_full="${processed}/${project}/${sample_name}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
	fi
	#echo "${ARDB_full}"
	while IFS= read -r line; do
		# exit if no genes were found for the sample
		if [[ -z "${line}" ]] || [[ "${line}" == *"No anti-microbial genes were found"* ]]; then
			break
		fi
		IFS='	' read -r -a ar_line <<< "$line"
		#echo "ls-${line_squeezed} in ${sample_name}.ResGANNOT.${2}_${3}_sstar_summary.txt"
		# gets the match length of the gene to the database
		#length_1=$(echo ${line} | cut -d$'\t' -f8)
		# gets the total length of the gene from the database
		#length_2=$(echo ${line} | cut -d$'\t' -f9)
		length_1="${ar_line[7]}"
		length_2="${ar_line[8]}"
		percent_ID="${ar_line[6]}"
		percent_length="${ar_line[9]}"
		conferred=$(echo "${ar_line[1]}" | cut -d'_' -f1)
		#if [[ "${ar_line[3]}" == *"trunc"* ]] || [[ "${ar_line[3]}" == "trunc"* ]] || [[ "${ar_line[3]}" == *"trunc" ]] || [[ "${ar_line[3]}" == "trunc" ]]; then
		#	gene="TRUNC-${ar_line[4]}"
		#else
		gene="${ar_line[4]}"
		#fi
		#echo "norm:${sample_name}:${line}:${length_1}|${length_2}"
		#percent_ID=$(echo ${line} | cut -d$'\t' -f7)
		#percent_length=$(echo ${line} | cut -d$'\t' -f10)
		# gets the name of the gene
		#gene=$(echo ${line} | cut -d$'\t' -f5)
		#conferred=$(echo ${line} | cut -d$'\t' -f2 | cut -d'_' -f1)
		if [[ "${conferred}" == "macrolide_" ]] || [[ "${conferred}" == "macrolide," ]]; then
			conferred="macrolide_lincosamide_streptogramin_B"
		fi
		# gets the difference in length of the match vs total length
		#alength=$((length_2 - length_1))
#		echo "$sample_name --- l1=${length_1};l2=${length_2} --- ${gene}"
		#percent_length=$(( 100 * length_1 / length_2 ))
#		echo "cutoff length of ${length_2} is ${percent_length}"
		# Checks to see if it is one of the genes that have multiple representations. If found it is changed to match the one accepted version
		#echo "pf_l-${percent_length}>=${project_parser_Percent_length}---pf_id-${percent_ID}>=${project_parser_Percent_identity}"
		if [[ ${percent_length} -ge ${project_parser_Percent_length} ]] && [[ ${percent_ID} -ge ${project_parser_Percent_identity} ]] ; then
			if [[ -z "${oar_list}" ]]; then
			#	echo "First oar: ${gene}"
				oar_list="${gene}(${conferred})[${percent_ID}/${percent_length}]"
			else
				if [[ ${oar_list} == *"${gene}"* ]]; then
				#	echo "${gene} already found in ${oar_list}"
					:
				else
				#	echo "${gene} not found in ${oar_list}...adding it"
					oar_list="${oar_list},${gene}(${conferred})[${percent_ID}/${percent_length}]"
				fi
			fi
		# If length is less than predetermined minimum (90% right now) then the gene is added to a rejects list to show it was outside acceptable limits
		else
			echo -e "${project}\t${sample_name}\tfull_assembly\t${line}" >> ${output_directory}/${4}-csstar_rejects.txt
		fi
	done < ${ARDB_full}
	# Changes list names if empty
	if [[ -z "${oar_list}" ]]; then
		oar_list="No AR genes discovered"
	fi
	mlst=$(head -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst)
	mlst=$(echo "${mlst}" | cut -d'	' -f3)
	if [[ ! -f "${processed}/${project}/${sample_name}/${sample_name}.txt" ]]; then
		"${shareScript}/determine_taxID.sh" "${sample_name}" "${project}"
	fi
	tax_file="${processed}/${project}/${sample_name}/${sample_name}.txt"
	#echo "Looking at ${processed}/${project}/${sample_name}/${sample_name}.tax"
	genus=$(tail -2 "${processed}/${project}/${sample_name}/${sample_name}.tax"| head -n1 | cut -d'	' -f2)
	species=$(tail -1 "${processed}/${project}/${sample_name}/${sample_name}.tax" | cut -d'	' -f2)
	ANI="${genus} ${species}"
#	echo "${ANI}"
	echo -e "${project}\t${sample_name}\t${ANI}\t${mlst}\t${oar_list}" >> ${output_directory}/${4}-csstar_summary_full.txt

	#Adding in srst2 output internalSTOPcodon
	if [[ -s "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${resGANNOT_srst2_filename}_srst2__results.txt" ]]; then
		srst2_results=""
		while IFS= read -r line; do
			echo "Start"
			gene=$(echo "${line}" | cut -d'	' -f3)
			#ODD WAY to do this right now, must look into later, but
			confers=$(echo "${line}" | cut -d'	' -f14 | cut -d';' -f3)
			echo "${gene}-${confers}"
			if [[ "${confers}" = "annotation" ]]; then
				continue
			fi
			if [[ -z "${confers}" ]]; then
				if [[ ! -z ${gene} ]]; then
					if [[ "${gene,,}" == "agly_flqn" ]]; then
						confers="aminoglycoside_and_fluoroquinolone_resistance"
					elif [[ "${gene,,}" == "tetracenomycinc" ]]; then
						confers="tetracenomycinC_resistance"
					else
						confers=${groups[${gene:0:3}]}
					fi
				fi
			fi
			confers=${confers//_resistance/}
			allele=$(echo "${line}" | cut -d'	' -f4 | cut -d'_' -f1)
			if [[ "${allele}" = "Zn-dependent" ]]; then
				allele="${allele}_hydrolase"
			fi
			coverage=$(echo "${line}" | cut -d'	' -f5)
			depth=$(echo "${line}" | cut -d'	' -f6)
			diffs=$(echo "${line}" | cut -d'	' -f7)
			if [[ ${diffs} == *"trunc"* ]]; then
				allele="TRUNC-${allele}"
			fi
			uncertainty=$(echo "${line}" | cut -d'	' -f8)
			divergence=$(echo "${line}" | cut -d'	' -f9)
			``
			length=$(echo "${line}" | cut -d'	' -f10)
			percent_length=$(echo "$coverage / 1" | bc)
			if [[ "${divergence}" = "0.0" ]]; then
				percent_ID=100
			else
				percent_ID=$(echo "100 - (($divergence + 1) / 1)" | bc)
			fi
			echo "${allele}/${coverage}/${depth}/${diffs}/${uncertainty}/${divergence}/${length}/${percent_ID}/${percent_length}"
			if [[ "${percent_ID}" -gt 95 ]] && [[ "${percent_length}" -gt 90 ]]; then
				info_line="${allele}(${confers})[${percent_ID}/${percent_length}]"
				if [[ -z "${srst2_results}" ]]; then
					srst2_results=${info_line,,}
				else
					srst2_results="${srst2_results},${info_line,,}"
				fi
			else
				if [[ ${line} = "Sample	DB	gene"* ]]; then
					:
				else
					echo ${line} >> ${output_directory}/${4}-srst2_rejects.txt
				fi
			fi
		done < "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${resGANNOT_srst2_filename}_srst2__results.txt"
		echo "Test1"
		if [[ -z "${srst2_results}" ]]; then
			echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${4}-srst2.txt
		else
			echo "${project}	${sample_name}	${srst2_results}" >> ${output_directory}/${4}-srst2.txt
		fi
	else
		echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${4}-srst2.txt
	fi

#Test
echo "Test"


	if [[ "${has_plasmidAssembly}" = "true" ]]; then
		# Repeat the c-sstar output organization of the plasmidAssembly
		oar_list=""
		# Looks at all the genes found on the plasmid assembly for a sample
		if [[ -f "${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt" ]]; then
			ARDB_plasmid="${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
		else
			${shareScript}/run_c-sstar_on_single.sh "${sample_name}" "${gapping}" "${sim}" "${project}" "--plasmid"
			ARDB_plasmid="${processed}/${project}/${sample_name}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
		fi
		while IFS= read -r line; do
			# exit if no genes were found for the sample
			if [[ "${line}" == *"No anti-microbial genes were found"* ]]; then
				break
			fi
			IFS='	' read -r -a ar_line <<< "$line"
			# gets the match length of the gene to the database
			#length_1=$(echo ${line} | cut -d'	' -f8)
			# gets the total length of the gene from the database
			#length_2=$(echo ${line} | cut -d'	' -f9)
	#		echo "plas:${sample_name}:${line_squeezed}:${length_1}|${length_2}"
			#percent_ID=$(echo ${line} | cut -d'	' -f7)
			#percent_length=$(echo ${line} | cut -d'	' -f10)
			# gets the name of the gene
			#gene=$(echo ${line} | cut -d'	' -f5)
			#conferred=$(echo ${line} | cut -d'	' -f2 | cut -d'_' -f1)
			length_1="${ar_line[7]}"
			length_2="${ar_line[8]}"
			percent_ID="${ar_line[6]}"
			percent_length="${ar_line[9]}"
			conferred=$(echo "${ar_line[1]}" | cut -d'_' -f1)
			gene="pla-${ar_line[4]}"
			if [[ "${ar_line[3]}" == *"trunc"* ]] || [[ "${ar_line[3]}" == "trunc"* ]] || [[ "${ar_line[3]}" == *"trunc" ]] || [[ "${ar_line[3]}" == "trunc" ]]; then
				gene="TRUNC-pla-${ar_line[4]}"
			else
				gene="pla-${ar_line[4]}"
			fi
			if [[ "${conferred}" == "macrolide," ]]; then
				conferred="macrolide, lincosamide, streptogramin_B"
			fi
			# gets the difference in length of the match vs total length
			#alength=$((length_2 - length_1))
	#		echo "$sample_name --- l1=${length_1};l2=${length_2} --- ${gene}"
			#percent_length=$(( 100 * length_1 / length_2 ))
	#		echo "cutoff length of ${length_2} is ${percent_length}"
			# Checks to see if it is one of the genes that have multiple representations. If found it is changed to match the one accepted version
			if [[ ${percent_length} -ge ${project_parser_Percent_length} ]] && [[ ${percent_ID} -ge ${project_parser_plasmid_Percent_identity} ]] ; then
				if [[ -z "${oar_list}" ]]; then
				#	echo "First oar: ${gene}"+
					oar_list="${gene}(${conferred})[${percent_ID}/${percent_length}]"
				else
					if [[ ${oar_list} == *"${gene}"* ]]; then
					#	echo "${gene} already found in ${oar_list}"
						:
					else
					#	echo "${gene} not found in ${oar_list}...adding it"
						oar_list="${oar_list},${gene}(${conferred})[${percent_ID}/${percent_length}]"
					fi
				fi
			# If length is less than predetermined minimum (90% right now) then the gene is added to a rejects list to show it was outside acceptable limits
			else
				echo -e "${project}\t${sample_name}\t${line}" >> ${output_directory}/${4}-csstar_rejects_plasmids.txt
			fi
		done < ${ARDB_plasmid}
		if [[ -z "${oar_list}" ]]; then
			oar_list="No AR genes discovered"
		fi
		echo -e "${project}\t${sample_name}\t${oxa_list}\t${oar_list}" >> ${output_directory}/${4}-csstar_summary_plasmid.txt
	fi


	# Goes through the plasmid file of the sample and adds all found plasmids to the summary file
	#echo "Starting plasmid extraction"
	if [[ -f ${processed}/${project}/${sample_name}/plasmid/${sample_name}_results_table_summary.txt ]]; then
		#echo "Found plasmid file"
		:
	fi
	full_contigs=">"
	full_contigs=$(grep -c ${full_contigs} "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta")
	added=0
	while IFS= read -r plasmid; do
		line_in=$(echo ${plasmid} | cut -d' ' -f1)
		if [[ "${line_in}" = "No" ]] || [[ "${line_in}" = "Enterococcus,Streptococcus,Staphylococcus" ]] || [[ "${line_in}" = "Enterobacteriaceae" ]] || [[ "${line_in}" = "Plasmid" ]]; then
			# echo "Not using line: $plasmid"
			:
		else
			echo -e "${project}\t${sample_name}\tfull_assembly\t${plasmid}" >> ${output_directory}/${4}-plasmid_summary.txt
			added=1
		fi
	done < ${processed}/${project}/${sample_name}/plasmid/${sample_name}_results_table_summary.txt
	if [[ "${added}" -eq 0 ]]; then
		echo -e "${project}\t${sample_name}\tfull_assembly\tNo_Plasmids_Found\t${full_contigs}_contigs-${components}_components" >> ${output_directory}/${4}-plasmid_summary.txt
	fi
	plas_contigs=">"
	plas_contigs=$(grep -c ${plas_contigs} "${OUTDATADIR}/plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta")
	components=-1
	while IFS= read -r contigs; do
		if [[ "${contigs}" = ">"* ]]; then
			components_temp=$(echo "${contigs}" | cut -d'_' -f8)
			if [[ ${components_temp} -gt ${components} ]]; then
				components="${components_temp}"
			fi
		fi
	done < ${processed}/${project}/${sample_name}/plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta
	components=$(( components + 1 ))
	while IFS= read -r plasmid; do
		line_in=$(echo ${plasmid} | cut -d' ' -f1)
		if [[ "${line_in}" = "No" ]] || [[ "${line_in}" = "Enterococcus,Streptococcus,Staphylococcus" ]] || [[ "${line_in}" = "Enterobacteriaceae" ]] || [[ "${line_in}" = "Plasmid" ]]; then
			:
		else
			echo -e "${project}\t${sample_name}\tplasmid_assembly\t${plasmid}" >> ${output_directory}/${4}-plasmid_summary.txt
			added=1
		fi
	done < ${processed}/${project}/${sample_name}/plasmid_on_plasmidAssembly/${sample_name}_results_table_summary.txt

	if [[ "${added}" -eq 0 ]]; then
		echo -e "${project}\t${sample_name}\tplasmid_assembly\tNo_Plasmids_Found\t${plas_contigs}_contigs-${components}_components" >> ${output_directory}/${4}-plasmid_summary.txt
	fi

	# Pulls MLST type for sample and adds it to the summary file
	if [[ -f "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" ]]; then
		mlst=$(head -n 1 ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst)
		mlst=$(echo "${mlst}" | cut -d'	' -f3-)
	else
		mlst="N/A"
	fi
	echo -e "${project}\t${sample_name}\t${mlst}" >> ${output_directory}/${4}-mlst_summary.txt
	#echo -e "${sample_name}\t${mlst}" >> ${share}/${4}-mlst_summary.txt

done < ${1}

python3 "${shareScript}/project_parser.py" "${output_directory}/${4}-csstar_summary_full.txt" "${output_directory}/${4}-plasmid_summary.txt" "${output_directory}/${4}_AR_plasmid_report.csv" "${output_directory}/${4}-srst2.txt"
