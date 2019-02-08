#!/bin/sh -l

#$ -o ablmq-cs.out
#$ -e ablmq-cs.err
#$ -N ablmq-cs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./abl_mass_qsub_csstar.sh path_to_list max_concurrent_submissions
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to ./abl_mass_qsub_csstar.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_csstar.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions output_directory_for_scripts"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"


# Create direcory to hold all temporary qsub scripts
counter=0
max_subs=${2}

# Set script directory
main_dir="${3}/csstar_alt_subs"
if [[ ! -d "${3}/csstar_alt_subs" ]]; then
	mkdir "${3}/csstar_alt_subs"
	mkdir "${3}/csstar_alt_subs/complete"
elif [[ ! -d  "${3}/csstar_alt_subs/complete" ]]; then
	mkdir "${3}/csstar_alt_subs/complete"
fi

# format name being extracted from alt database
main_dir="${3}/csstar_subs"
if [[ ! -d "${3}/csstar_subs" ]]; then
	mkdir -p "${3}/csstar_subs/complete"
elif [[ ! -d  "${3}/csstar_subs/complete" ]]; then
	mkdir "${3}/csstar_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit scripts to run default csstar on all samples on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	#rm -r ${processed}/${project}/${sample}/c-sstar/*20181003*
	#rm -r ${processed}/${project}/${sample}/c-sstar_plasmid/*20181003*
	echo ${counter}-${project}-${sample}
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		echo "Test"
		if [[ ${counter} -lt ${max_subs} ]]; then
			#if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				# Defaulting to gapped/98, change if you want to include user preferences
				echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g h \"${project}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				cd "${main_dir}"
				if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
				else
					if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
						qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
					else
						qsub -sync y "${main_dir}/csstn_${sample}_${start_time}.sh"
					fi
				fi
			#else
			#	echo "${project}/${sample} already has ${resGANNOT_srst2_filename}"
			#	echo "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
			#fi

			if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
				#if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${resGANNOT_srst2_filename}.gapped_40_sstar_summary.txt" ]]; then
					echo "Index below max submissions, submitting plasmid"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					# Defaulting to gapped/98, change if you want to include user preferences
					#echo -e "\"${shareScript}/fix_node_rename.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g o \"${project}\" \"--plasmid\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
					cd "${main_dir}"
					if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
					else
						qsub -sync y "${main_dir}/csstp_${sample}_${start_time}.sh"
					fi
				#else
				#	echo "${project}/${sample} already has ${resGANNOT_srst2_filename} PLASMID"
				#	echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
				#fi
			fi


		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			while :
			do
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
					#if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${resGANNOT_srst2_filename}.gapped_98_sstar_summary.txt" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						# Defaulting to gapped/98, change if you want to include user preferences
						echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g h \"${project}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
						else
							if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
								qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
							else
								qsub -sync y "${main_dir}/csstn_${sample}_${start_time}.sh"
							fi
						fi
					#else
					#	echo "${project}/${sample} already has ${resGANNOT_srst2_filename}"
					#	echo "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
					#fi

					if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
						#if [[ ! -f "${processed}/${project}/${sample}/c-sstar_plasmid/${sample}.${resGANNOT_srst2_filename}.gapped_40_sstar_summary.txt" ]]; then
							echo -e "#!/bin/bash -l\n" > "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -o csstp_${sample}.out" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -e csstp_${sample}.err" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -N csstp_${sample}"   >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -cwd"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "#$ -q short.q\n"  >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							# Defaulting to gapped/98, change if you want to include user preferences
							echo -e "\"${shareScript}/run_c-sstar_on_single.sh\" \"${sample}\" g o \"${project}\" \"--plasmid\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarp_complete.txt\"" >> "${main_dir}/csstp_${sample}_${start_time}.sh"
							cd ${main_dir}
							if [[ "${counter}" -lt "${last_index}" ]]; then
								qsub "${main_dir}/csstp_${sample}_${start_time}.sh"
							else
								qsub -sync y "${main_dir}/csstp_${sample}_${start_time}.sh"
							fi
						#else
						#	echo "${project}/${sample} already has ${resGANNOT_srst2_filename} PLASMID"
						#	echo "$(date)" > "${main_dir}/complete/${sample}_csstarp_complete.txt"
						#fi
					fi
					break
				else
					timer=$(( timer + 5 ))
					echo "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt not ready, sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of last sample
finish_counter=0
waiting_for_index=${last_index}
waiting_sample=$(echo "${arr[${last_index}]}" | cut -d'/' -f2)
timer=0
while :
do
	if [[ ${timer} -gt 3600 ]]; then
		echo "Timer exceeded limit of 3600 seconds = 60 minutes"
		break
	fi
	if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "All isolates completed"
		printf "%s %s" "Act_by_list.sh has completed ${2}" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
		exit 0
	else
		counter+$(( counter + 1))
		timer=$(( timer + 5 ))
		echo "sleeping for 5 seconds, so far slept for ${timer}"
		sleep 5
	fi
done
