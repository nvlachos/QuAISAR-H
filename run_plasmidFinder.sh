#!/bin/sh -l

#$ -o get_run_plasmidFinder.out
#$ -e get_run_plasmidFinder.err
#$ -N get_run_plasmidFinder
#$ -cwd
#$ -q all.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
. ${mod_changers}/pipeline_mods
. ${mod_changers}/list_modules.sh

#
# Will attempt to find any plasmids in sample
#
# Usage ./get_run_plasmidFinder.sh sample_name run_id force
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_plasmidFinder.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_plasmidFinder.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_plasmidFinder.sh  sample_name run_id output_folder(either plasmid or (DONT USE-plasmid_on_plasmidAssembly) (-i number_minimum_identity, optional) (-f to force against all databases, optional)"
	echo "Output by default is ${processed}/miseq_run_id/sample_name/plasmid"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty miseq_run_id supplied to run_plasmidFinder.sh, exiting"
	exit 1
elif [[ "${4}" == "-f" ]] || [[ "${6}" == "-f" ]]; then
	force="true"
fi

# Create output directory
if [[ ! -d ${processed}/${2}/${1}/${3} ]]; then
	echo "Making ${processed}/${2}/${1}/${3}"
	mkdir ${processed}/${2}/${1}/${3}
fi

# Set output directory
OUTDATADIR=${processed}/${2}/${1}/${3}
# Get proper input file based on output directory (whether it is full assembly or plasmid)
# if [[ "${3}" == "plasmid_on_plasmidAssembly" ]]; then
# 	outpath="plasmidAssembly/${1}_plasmid_scaffolds_trimmed.fasta"
# el
if [[ "${3}" == "plasmid" ]]; then
	outpath="Assembly/${1}_scaffolds_trimmed.fasta"
else
	echo "Non standard output location, using full assembly to find plasmids"
	outpath="Assembly/${1}_scaffolds_trimmed.fasta"
fi

# If flag was set to change % ID threshold then use it, otherwise use default from config.sh
if [[ "${4}" = "-i" ]]; then
	pl_id="${5}"
else
	pl_id="${plasmid_identity}"
fi


#If force flag is set, then run it against all databases
if [[ "${force}" == "true" ]]; then
	echo "Checking against ALL plasmids, but unlikely to find anything"
	plasmidfinder -i ${processed}/${2}/${1}/${outpath} -o ${OUTDATADIR} -k ${plasmidFinder_identity} -p enterobacteriaceae
	# Rename all files to include ID
	mv ${OUTDATADIR}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${1}_Hit_in_genome_seq_entero.fsa
	mv ${OUTDATADIR}/Plasmid_seq.fsa ${OUTDATADIR}/${1}_Plasmid_seq_enetero.fsa
	mv ${OUTDATADIR}/results.txt ${OUTDATADIR}/${1}_results_entero.txt
	mv ${OUTDATADIR}/results_tab.txt ${OUTDATADIR}/${1}_results_tab_entero.txt
	mv ${OUTDATADIR}/results_table.txt ${OUTDATADIR}/${1}_results_table_entero.txt
	plasmidfinder -i ${processed}/${2}/${1}/${outpath} -o ${OUTDATADIR} -k ${plasmidFinder_identity} -p gram_positive
	# Rename all files to include ID
	mv ${OUTDATADIR}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${1}_Hit_in_genome_seq_gramp.fsa
	mv ${OUTDATADIR}/Plasmid_seq.fsa ${OUTDATADIR}/${1}_Plasmid_seq_gramp.fsa
	mv ${OUTDATADIR}/results.txt ${OUTDATADIR}/${1}_results_gramp.txt
	mv ${OUTDATADIR}/results_tab.txt ${OUTDATADIR}/${1}_results_tab_gramp.txt
	mv ${OUTDATADIR}/results_table.txt ${OUTDATADIR}/${1}_results_table_gramp.txt
	cat	${OUTDATADIR}/${1}_results_table_gramp.txt ${OUTDATADIR}/${1}_results_table_entero.txt > ${OUTDATADIR}/${1}_results_table_summary.txt
# Else, if the force flag is not set, then TRY to limit search to family (it will still check against all if it does not match the family)
else
	# Checks to see if a post assembly kraken file is available to extract the family of the sample
	if [[ -f ${processed}/${2}/${1}/kraken/postAssembly/${1}_kraken_summary_assembled_BP_data.txt ]]; then
		#Extracts the 6th line from the kraken file containing all family information
		family=$(sed -n '6p' < ${processed}/${2}/${1}/kraken/postAssembly/${1}_kraken_summary_assembled.txt)
		genus=$(sed -n '7p' < ${processed}/${2}/${1}/kraken/postAssembly/${1}_kraken_summary_assembled.txt)
		#Extracts family name from line
		family=$(echo ${family}  | cut -d' ' -f4)
		genus=$(echo ${genus}  | cut -d' ' -f4)
		echo "${family}-${genus}"
		# If family is enterobacteriaceae, then run against that DB
		if [[ "${family,}" == "enterobacteriaceae" ]]; then
			echo "Checking against Enterobacteriaceae plasmids"
			plasmidfinder -i ${processed}/${2}/${1}/${outpath} -o ${OUTDATADIR} -k ${plasmidFinder_identity} -p enterobacteriaceae
			# Rename all files to include ID
			mv ${OUTDATADIR}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${1}_Hit_in_genome_seq_entero.fsa
			mv ${OUTDATADIR}/Plasmid_seq.fsa ${OUTDATADIR}/${1}_Plasmid_seq_enetero.fsa
			mv ${OUTDATADIR}/results.txt ${OUTDATADIR}/${1}_results_entero.txt
			mv ${OUTDATADIR}/results_tab.txt ${OUTDATADIR}/${1}_results_tab_entero.txt
			mv ${OUTDATADIR}/results_table.txt ${OUTDATADIR}/${1}_results_table_summary.txt
		# If family is staph, strp, or enterococcus, then run against the gram positive database
		elif [[ "${genus,}" == "staphylococcus" ]] || [[ "${3,}" == "streptococcus" ]] || [[ "${3,}" == "enterococcus" ]]; then
			echo "Checking against Staph, Strep, and Enterococcus plasmids"
			plasmidfinder -i ${processed}/${2}/${1}/${outpath} -o ${OUTDATADIR} -k ${plasmidFinder_identity} -p gram_positive
			# Rename all files to include ID
			mv ${OUTDATADIR}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${1}_Hit_in_genome_seq_gramp.fsa
			mv ${OUTDATADIR}/Plasmid_seq.fsa ${OUTDATADIR}/${1}_Plasmid_seq_gramp.fsa
			mv ${OUTDATADIR}/results.txt ${OUTDATADIR}/${1}_results_gramp.txt
			mv ${OUTDATADIR}/results_tab.txt ${OUTDATADIR}/${1}_results_tab_gramp.txt
			mv ${OUTDATADIR}/results_table.txt ${OUTDATADIR}/${1}_results_table_summary.txt
		# Family is not one that has been designated by the creators of plasmidFinder to work well, but still attempting to run against both databases
		else
			echo "Checking against ALL plasmids, but unlikely to find anything"
			plasmidfinder -i ${processed}/${2}/${1}/${outpath} -o ${OUTDATADIR} -k ${plasmidFinder_identity} -p enterobacteriaceae
			# Rename all files to include ID
			mv ${OUTDATADIR}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${1}_Hit_in_genome_seq_entero.fsa
			mv ${OUTDATADIR}/Plasmid_seq.fsa ${OUTDATADIR}/${1}_Plasmid_seq_enetero.fsa
			mv ${OUTDATADIR}/results.txt ${OUTDATADIR}/${1}_results_entero.txt
			mv ${OUTDATADIR}/results_tab.txt ${OUTDATADIR}/${1}_results_tab_entero.txt
			mv ${OUTDATADIR}/results_table.txt ${OUTDATADIR}/${1}_results_table_entero.txt
			plasmidfinder -i ${processed}/${2}/${1}/${outpath} -o ${OUTDATADIR} -k ${plasmidFinder_identity} -p gram_positive
			# Rename all files to include ID
			mv ${OUTDATADIR}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${1}_Hit_in_genome_seq_gramp.fsa
			mv ${OUTDATADIR}/Plasmid_seq.fsa ${OUTDATADIR}/${1}_Plasmid_seq_gramp.fsa
			mv ${OUTDATADIR}/results.txt ${OUTDATADIR}/${1}_results_gramp.txt
			mv ${OUTDATADIR}/results_tab.txt ${OUTDATADIR}/${1}_results_tab_gramp.txt
			mv ${OUTDATADIR}/results_table.txt ${OUTDATADIR}/${1}_results_table_gramp.txt
			cat	${OUTDATADIR}/${1}_results_table_gramp.txt ${OUTDATADIR}/${1}_results_table_entero.txt > ${OUTDATADIR}/${1}_results_table_summary.txt

		fi
	# No assembly file exists and cannot be used to determine family of sample
	else
		echo "Cant guess the genus of the sample, please try again with the force option or check the contents of the postAssembly kraken folder for a proper summary file (${processed}/${2}/${1}/kraken/postAssembly/${1}_kraken_summary_assembled_BP_data.txt)"
	fi
fi
