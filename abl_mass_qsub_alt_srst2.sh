#!/bin/sh -l

#$ -o ablmq_alts.out
#$ -e ablmq_alts.err
#$ -N ablmq_alts
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"


#
# Usage ./abl_mass_qsub_alt_srst2.sh path_to_list max_concurrent_submissions path_to_alt_database output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_alt_srst2.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions path_to_alt_database output_directory_for_scripts"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ ! -f "${3}" ]] || [[ -z "${3}" ]]; then
	echo "${3} (alt_db) does not exist...exiting"
	exit 1
elif [[ -z  "${4}" ]]; then
	echo "outdir  parameter for scripts is empty...exiting"
	exit 1
elif [[ -z "${5}" ]]; then
	echo "clobberness is empty...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${5}" != "keep" ]] && [[ "${5}" != "clobber" ]]; then
	echo "Clobberness was not input correctly [keep|clobber]...exiting"
	exit 1
else
	clobberness="${5}"
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [ "$line" ];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

# Create counter and set max number of concurrent submissions
counter=0
max_subs=${2}

# Set script directory
main_dir="${4}/srst2_alt_subs"
if [[ ! -d "${4}/srst2_alt_subs" ]]; then
	mkdir "${4}/srst2_alt_subs"
	mkdir "${4}/srst2_alt_subs/complete"
elif [[ ! -d  "${4}/srst2_alt_subs/complete" ]]; then
	mkdir "${4}/srst2_alt_subs/complete"
fi

# Format name of DB used from filename for srst2
alt_DB_path=${3}
alt_DB=$(echo ${alt_DB_path##*/} | cut -d'.' -f1)
alt_DB=${alt_DB//_srst2/}
echo ${alt_DB}
start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit srst2 qsub scripts for each isolate in the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	# Fix any improperly named (older style) files
	if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__genes__${alt_DB}_srst2__results.txt" ]]; then
			rm "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__genes__${alt_DB}_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt"
			if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
				rm "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${alt_DB}_srst2__results.txt"
			fi
			continue
	elif [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
		rm "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt"
		continue
	fi
	# Remove old data if clobber is set
	if [[ "${clobberness}" == "clobber" ]]; then
		if [[ -f "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__genes__${alt_DB}_srst2__results.txt" ]]; then
			rm "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__genes__${alt_DB}_srst2__results.txt"
		fi
		if [[ "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
			rm "${processed}/${project}/${sample_name}/srst2/${sample_name}_${alt_DB}__fullgenes__${alt_DB}_srst2__results.txt"
		fi
	fi
	echo ${counter}
	# If counter is below max number of concurrent submissions...submit away
	if [ ${counter} -lt ${max_subs} ]; then
		# Checks if old results exist, skips if so
		if [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt" ]] || [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
			echo -e "cd ${shareScript}" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "\"${shareScript}/run_srst2AR_altDB.sh\" \"${sample}\" \"${project}\" \"${alt_DB_path}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2AR_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			cd "${main_dir}"
			qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
		# Old data exists, skipping
		else
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2AR_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo "${project}/${sample} already has ${alt_DB}"
		fi
	# Counter is above max xubmissions, must wait for "slot" to open up
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
		# Loop to wait until slot opens up
		while :
		do
			# if max time limit has been reached, then exit
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			# Check if "waiting" sample has completed analysis
			if [ -f "${main_dir}/complete/${waiting_sample}_srst2AR_complete.txt" ]; then
				# Check for old data, skip if present
				if [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${alt_DB}_srst2__results.txt" ]] || [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${alt_DB}_srst2__results.txt" ]]; then
					echo "${waiting_sample} has completed, starting ${sample}"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
					echo -e "cd ${shareScript}" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/run_srst2AR_altDB.sh\" \"${sample}\" \"${project}\" \"${alt_DB_path}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2AR_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					cd "${main_dir}"
					qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
				# Old data exists, skipping
				else
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2AR_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo "${project}/${sample} already has 20180608"
				fi
				break
			# Wait 5 seconds before checking if "waiting" sample is finished yet
			else
				timer=$(( timer + 5 ))
				echo "sleeping for 5 seconds, so far slept for ${timer}"
				sleep 5
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_srst2AR_complete.txt" ]]; then
		echo "${item} is complete"
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_srst2AR_complete.txt" ]]; then
					echo "${item} is complete"
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
printf "%s %s" "abl_mass_qsub_alt_srst2.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_alt_srst2.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
