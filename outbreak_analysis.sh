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
# Usage ./outbreak_analysis.sh path_to_list gapped/ungapped (analysis ran) identity (80/95/98/99/100) analysis_identifier(e.g. outbreak identifier) output_directory(will create a folder at this location with nam of analysis_identifier) plasmid_identity(optional)
#
# Pulls out MLST, AR genes, and plasmid repicons and creates a mashtree for the listed samples and consolidates them into one sheet
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_csstar_proj_parser.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./outbreak_analysis.sh path_to_list_file gapped/ungapped 80/95/98/99/100 output_prefix output_directory plasmid_identity_cutoff(optional, default = 40)"
	exit 0
elif [[ ! -f ${1} ]]; then
	echo "list does not exit...exiting"
	exit 1
fi
# Checks that the gapping is set to one of the csstar presets
if [[ "${2}" != "gapped" ]] && [[ "${2}" != "ungapped" ]]; then
	echo "gapping does not equal; gapped or ungapped...exiting"
	exit 1
fi
# Checks that value given for % Identity is one of the presets for csstar
if [[ "${3}" != 80 ]] && [[ "${3}" != 95 ]] && [[ "${3}" != 98 ]] && [[ "${3}" != 99 ]] && [[ "${3}" != 100 ]]; then
	echo "Identity is not one of the presets for csstar and therefore will fail, exiting..."
	exit 1
fi
# Creates the output directory if it does not exist
output_directory=${5}/${4}
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

# Creates a dictionary to match genes to AR conferred when parsing srst files
declare -A groups
echo ""
echo "Creating AR lookup list from ${local_DBs}/star/group_defs.txt"
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
done < "${local_DBs}/star/group_defs.txt"

# Set defaults for checking if all isolates have been compared to the newest ResGANNOT DB file . If any isolates are not up-to-date, they will be submitted with the appropriate abl_mass_qsub.
run_csstar="false"
run_srst2="false"
> "${output_directory}/${4}_csstar_todo.txt"
> "${output_directory}/${4}_srst2_todo.txt"

# Check that each isolate has been comparde to the newest ResGANNOT DB file
echo -e "\nMaking sure all isolates use the latest AR Database - ${resGANNOT_srst2_filename}\n"
while IFS= read -r line; do
	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	OUTDATADIR="${processed}/${project}/${sample_name}"
	#echo "checking for ${OUTDATADIR}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
	if [[ -s "${OUTDATADIR}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt" ]];
	then
		echo "${project}/${sample_name} has newest ResGANNOT for normal csstar already"
	else
		echo "${project}/${sample_name}" >> "${output_directory}/${4}_csstar_todo.txt"
		run_csstar="true"
	fi
	#echo "checking for ${OUTDATADIR}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
	if [[ -s "${OUTDATADIR}/plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta" ]]; then
		if [[ -s "${OUTDATADIR}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt" ]] || [[ -s "${OUTDATADIR}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_40_sstar_summary.txt" ]]; then
			echo "${project}/${sample_name} has newest ResGANNOT for plasmid csstar already"
		else
			echo "${project}/${sample_name}" >> "${output_directory}/${4}_csstar_todo.txt"
			sort -u "${output_directory}/${4}_csstar_todo.txt" > "${output_directory}/${4}_csstar_todo_no_dups.txt"
			cp "${output_directory}/${4}_csstar_todo_no_dups.txt" "${output_directory}/${4}_csstar_todo.txt"
			run_csstar="true"
		fi
	else
		echo "No plasmid Assembly found, no need for csstar plasmid"
	fi
	#echo "checking for ${OUTDATADIR}/srst2/${sample_name}__genes__${resGANNOT_srst2_filename}_srst2__results.txt"
	if [[ -s ${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq ]] && [[ -s ${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq ]] || [[ -s ${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz ]] && [[ -s ${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
		#echo "FASTQs exist"
		if [[ -f "${OUTDATADIR}/srst2/${sample_name}__fullgenes__${resGANNOT_srst2_filename}_srst2__results.txt" ]] || [[ -f "${OUTDATADIR}/srst2/${sample_name}__genes__${resGANNOT_srst2_filename}_srst2__results.txt" ]]; then
				echo "${project}/${sample_name} has newest ResGANNOT for srst2 already"
			else
				echo "${project}/${sample_name}" >> "${output_directory}/${4}_srst2_todo.txt"
				run_srst2="true"
		fi
	fi
done < ${1}

# Creating mashtree of all isolates in list
echo "Creating mashtree of all samples"
${shareScript}/mashtree_of_list.sh "${1}" "${output_directory}/mashtree" "${4}"
cp "${output_directory}/mashtree/${4}.dnd" "${output_directory}/${4}.nwk"
sed -i "s/_scaffolds_trimmed//g" ${thing}
rm -r ${output_directory}/mashtree

# Submits the list of isolates that need the newest ResGANNOT file for csstar
if [[ "${run_csstar}" = "true" ]]; then
	echo "Submitting list for csstar qsub analysis"
	qsub ./abl_mass_qsub_csstar.sh "${output_directory}/${4}_csstar_todo.txt" 25
fi
# Submits the list of isolates that need the newest ResGANNOT file for srst2
if [[ "${run_srst2}" = "true" ]]; then
	echo "Submitting list for srst2 qsub analysis"
	qsub -sync y ./abl_mass_qsub_srst2.sh "${output_directory}/${4}_srst2_todo.txt" 25
fi

# # Loop through and extracts and formats AR genes found in all isolates, as well as the primary MLST type and plasmid replicons. Each are output to separate files. Any AR genes that do not meet the length or % identity are copied to the rejects file.
while IFS= read -r line; do
 	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
 	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
 	OUTDATADIR="${processed}/${project}/${sample_name}"
	sample_index=0
	oar_list=""
	# Looks at all the genes found for a sample
	#echo "looking for ${OUTDATADIR}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
	if [[ -f "${OUTDATADIR}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt" ]]; then
		ARDB_full="${OUTDATADIR}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
	else
		echo "IT STILL thinks it needs to run ${sample_name} through normal csstar"
		#${shareScript}/run_c-sstar_on_single.sh "${sample_name}" "${gapping}" "${sim}" "${project}"
		#ARDB_full="${OUTDATADIR}/c-sstar/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
	fi
	#echo "${ARDB_full}"
	# Extracts all AR genes from normal csstar output file and creates a lits of all genes that pass the filtering steps
	while IFS= read -r line; do
		# exit if no genes were found for the sample
		if [[ -z "${line}" ]] || [[ "${line}" == *"No anti-microbial genes were found"* ]]; then
			break
		fi
		IFS='	' read -r -a ar_line <<< "$line"
		length_1="${ar_line[7]}"
		length_2="${ar_line[8]}"
		percent_ID="${ar_line[6]}"
		percent_length="${ar_line[9]}"
		conferred=$(echo "${ar_line[1]}" | rev | cut -d'_' -f2- | rev)
		gene="${ar_line[4]}"
		# Ensure that the gene passes % identity and % length threhsolds for reporting
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

	# Extracts the MLST type
	mlst=$(head -n1 ${OUTDATADIR}/MLST/${sample_name}.mlst)
	mlst=$(echo "${mlst}" | cut -d'	' -f3)

	# Extracts taxonomic info
	if [[ ! -f "${OUTDATADIR}/${sample_name}.txt" ]]; then
		"${shareScript}/determine_taxID.sh" "${sample_name}" "${project}"
	fi
	tax_file="${OUTDATADIR}/${sample_name}.txt"
	#echo "Looking at ${OUTDATADIR}/${sample_name}.tax"
	genus=$(tail -2 "${OUTDATADIR}/${sample_name}.tax"| head -n1 | cut -d'	' -f2)
	species=$(tail -1 "${OUTDATADIR}/${sample_name}.tax" | cut -d'	' -f2)
	ANI="${genus} ${species}"
#	echo "${ANI}"
# Print all extracted info to primary file
	echo -e "${project}\t${sample_name}\t${ANI}\t${mlst}\t${oar_list}" >> ${output_directory}/${4}-csstar_summary_full.txt

	# Adding in srst2 output in a similar fashion as to how the csstar genes are output to the file.
	if [[ -s "${OUTDATADIR}/srst2/${sample_name}__fullgenes__${resGANNOT_srst2_filename}_srst2__results.txt" ]]; then
		srst2_results=""
		while IFS= read -r line; do
		#	echo "Start"
			gene=$(echo "${line}" | cut -d'	' -f3)
			#ODD WAY to do this right now, must look into later, but
			confers=$(echo "${line}" | cut -d'	' -f14 | cut -d';' -f3)
		#	echo "${gene}-${confers}"
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
		#	echo "${allele}/${coverage}/${depth}/${diffs}/${uncertainty}/${divergence}/${length}/${percent_ID}/${percent_length}"
		# Filter genes based on thresholds for length and percent identity
			if [[ "${percent_ID}" -gt ${project_parser_Percent_identity} ]] && [[ "${percent_length}" -gt ${project_parser_Percent_length} ]]; then
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
		done < "${OUTDATADIR}/srst2/${sample_name}__fullgenes__${resGANNOT_srst2_filename}_srst2__results.txt"
		#echo "Test1"
		if [[ -z "${srst2_results}" ]]; then
			echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${4}-srst2.txt
		else
			echo "${project}	${sample_name}	${srst2_results}" >> ${output_directory}/${4}-srst2.txt
		fi
	else
		echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${4}-srst2.txt
	fi

#Test
#echo "Test"

	# Parse through plasmid Assembly, although it is not used in the final report
	if [[ "${has_plasmidAssembly}" = "true" ]]; then
		# Repeat the c-sstar output organization of the plasmidAssembly
		oar_list=""
		# Looks at all the genes found on the plasmid assembly for a sample
		if [[ -f "${OUTDATADIR}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt" ]]; then
			ARDB_plasmid="${OUTDATADIR}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
		else
			echo "It STILL STILL thinks it needs to put ${sample_name} trhough plasmid csstar"
			#${shareScript}/run_c-sstar_on_single.sh "${sample_name}" "${gapping}" "${sim}" "${project}" "--plasmid"
			#ARDB_plasmid="${OUTDATADIR}/c-sstar_plasmid/${sample_name}.${resGANNOT_srst2_filename}.${2}_${3}_sstar_summary.txt"
		fi
		while IFS= read -r line; do
			# exit if no genes were found for the sample
			if [[ "${line}" == *"No anti-microbial genes were found"* ]]; then
				break
			fi
			IFS='	' read -r -a ar_line <<< "$line"
			length_1="${ar_line[7]}"
			length_2="${ar_line[8]}"
			percent_ID="${ar_line[6]}"
			percent_length="${ar_line[9]}"
			conferred=$(echo "${ar_line[1]}" | cut -d'_' -f1)
			gene="pla-${ar_line[4]}"
			if [[ "${conferred}" == "macrolide," ]]; then
				conferred="macrolide, lincosamide, streptogramin_B"
			fi
			# Checks to see if gene passes the threshold rquirements for identity and length
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
		# Adds generic output saying nothing was found if the list was empty
		if [[ -z "${oar_list}" ]]; then

			oar_list="No AR genes discovered"
		fi
		# Adds info to plasmid csstar summary file
		echo -e "${project}\t${sample_name}\t${oxa_list}\t${oar_list}" >> ${output_directory}/${4}-csstar_summary_plasmid.txt
	fi


	# Goes through the plasmid file of the sample and adds all found plasmid replicons to the summary file
	#echo "Starting plasmid extraction"
	if [[ -f ${OUTDATADIR}/plasmid/${sample_name}_results_table_summary.txt ]]; then
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
	done < ${OUTDATADIR}/plasmid/${sample_name}_results_table_summary.txt
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
	done < ${OUTDATADIR}/plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta
	components=$(( components + 1 ))
	while IFS= read -r plasmid; do
		line_in=$(echo ${plasmid} | cut -d' ' -f1)
		if [[ "${line_in}" = "No" ]] || [[ "${line_in}" = "Enterococcus,Streptococcus,Staphylococcus" ]] || [[ "${line_in}" = "Enterobacteriaceae" ]] || [[ "${line_in}" = "Plasmid" ]]; then
			:
		else
			echo -e "${project}\t${sample_name}\tplasmid_assembly\t${plasmid}" >> ${output_directory}/${4}-plasmid_summary.txt
			added=1
		fi
	done < ${OUTDATADIR}/plasmid_on_plasmidAssembly/${sample_name}_results_table_summary.txt

	if [[ "${added}" -eq 0 ]]; then
		echo -e "${project}\t${sample_name}\tplasmid_assembly\tNo_Plasmids_Found\t${plas_contigs}_contigs-${components}_components" >> ${output_directory}/${4}-plasmid_summary.txt
	fi

	# Pulls MLST type for sample and adds it to the summary file
	if [[ -f "${OUTDATADIR}/MLST/${sample_name}.mlst" ]]; then
		mlst=$(head -n 1 ${OUTDATADIR}/MLST/${sample_name}.mlst)
		mlst=$(echo "${mlst}" | cut -d'	' -f3-)
	else
		mlst="N/A"
	fi
	echo -e "${project}\t${sample_name}\t${mlst}" >> ${output_directory}/${4}-mlst_summary.txt
	#echo -e "${sample_name}\t${mlst}" >> ${share}/${4}-mlst_summary.txt
done < ${1}

# Calls script that sorts and formats all isolates info into a atrix for easy viewing
python3 "${shareScript}/project_parser.py" "${output_directory}/${4}-csstar_summary_full.txt" "${output_directory}/${4}-plasmid_summary.txt" "${output_directory}/${4}_AR_plasmid_report.csv" "${output_directory}/${4}-srst2.txt"
