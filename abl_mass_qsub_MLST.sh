#!/bin/sh -l

#$ -o ablmq-mlst.out
#$ -e ablmq-mlst.err
#$ -N ablmq-mlst
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi

. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#List all currently loaded modules
#. ./module_changers/list_modules.sh

#
# Usage ./abl_mass_qsub_MLST.sh path_to_list max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to ./abl_mass_qsub_MLST.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_MLST.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber]"
	exit 1


# Check that clobberness is a valid option
if [[ "${4}" != "keep" ]] && [[ "${4}" != "clobber" ]]; then
	echo "Clobberness was not input correctly, be sure to add keep or clobber as 5th parameter...exiting"
	exit 1
else
	clobberness="${4}"
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
main_dir="${3}/mlst_subs"
if [[ ! -d "${3}/mlst_subs" ]]; then
	mkdir "${3}/mlst_subs"
	mkdir "${3}/mlst_subs/complete"
elif [[ ! -d "${3}/mlst_subs/complete" ]]; then
	mkdir "${3}/mlst_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub script to determine MlST of all samples on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" == "clobber" ]]; then
		rm "${processed}/${project}/${sample}/${sample}.mlst"
	fi
	# Check if there is an acceptable assembly file to use with MLST
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		# Check if counter is below max sub limit
		if [[ ${counter} -lt ${max_subs} ]]; then
			# Check for old data, skip if present
			if [[ ! -f "${processed}/${project}/${sample}/MLST/${sample}.mlst" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/mlst_${sample}_${start_time}.sh"
				echo -e "#$ -o mlst_${sample}.out" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
				echo -e "#$ -e mlst_${sample}.err" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
				echo -e "#$ -N mlst_${sample}"   >> "${main_dir}/mlst_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/mlst_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/mlst_${sample}_${start_time}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_MLST.sh\" \"${sample}\" \"${project}\" \"-f\" \"abaumannii\"" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_mlst_complete.txt\"" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
				cd "${main_dir}"
				#if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/mlst_${sample}_${start_time}.sh"
				#else
				#	qsub -sync y "${main_dir}/mlst_${sample}_${start_time}.sh"
				#fi
			else
				echo "${project}/${sample} already has mlst"
				echo "$(date)" > "${main_dir}/complete/${sample}_mlst_complete.txt"
			fi
		# Counter is over limit, must wait until slot becomes available
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			# Loop to check if "waiting" sample is complete yet
			while :
			do
				# Check if timer is above max time allowed
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check if old data exists, skip if so
				if [[ -f "${main_dir}/complete/${waiting_sample}_mlst_complete.txt" ]]; then
					if [[ ! -f "${processed}/${project}/${sample}/MLST/${sample}.mlst" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/mlst_${sample}_${start_time}.sh"
						echo -e "#$ -o mlst_${sample}.out" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
						echo -e "#$ -e mlst_${sample}.err" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
						echo -e "#$ -N mlst_${sample}"   >> "${main_dir}/mlst_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/mlst_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/mlst_${sample}_${start_time}.sh"
						echo -e "cd ${shareScript}" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_MLST.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_mlst_complete.txt\"" >> "${main_dir}/mlst_${sample}_${start_time}.sh"
						cd "${main_dir}"
						#if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/mlst_${sample}_${start_time}.sh"
						#else
						#	qsub -sync y "${main_dir}/mlst_${sample}_${start_time}.sh"
						#fi
					# Old data exists skipping
					else
						echo "${project}/${sample} already has mlst summary"
						echo "$(date)" > "${main_dir}/complete/${sample}_mlst_complete.txt"
					fi
					break
				# Wait 5 seconds and then check if "waiting" sample completed yet
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	# No assembly file exists to run mlst on
	else
		echo "${project}/${sample} does not have assembly file to run mlst on"
		echo "$(date)" > "${main_dir}/complete/${sample}_mlst_complete.txt"
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_mlst_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/mlst_${sample}.out" ]]; then
			mv "${shareScript}/mlst_${sample}.out" ${main_dir}
		fi
		if [[ -f "${shareScript}/mlst_${sample}.err" ]]; then
 			mv "${shareScript}/mlst_${sample}.err" ${main_dir}
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_mlst_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/mlst_${sample}.out" ]]; then
						mv "${shareScript}/mlst_${sample}.out" ${main_dir}
					fi
					if [[ -f "${shareScript}/mlst_${sample}.err" ]]; then
			 			mv "${shareScript}/mlst_${sample}.err" ${main_dir}
					fi
					break
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
		done
	fi
done

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "abl_mass_qsub_mlst.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_mlst.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
