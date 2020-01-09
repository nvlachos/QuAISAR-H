#!/bin/sh -l

#$ -o ablmq_GAMA.out
#$ -e ablmq_GAMA.err
#$ -N ablmq_GAMA
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
# Usage ./abl_mass_qsub_GAMA.sh path_to_list max_concurrent_submissions path_to_alternate_DB output_directory_for_scripts clobberness[keep|clobber]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_GAMA.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions path_to_alt_database output_directory_for_scripts clobberness[keep|clobber]"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${4}" ]]; then
	echo "Output directory parameter is empty...exiting"
	exit 1
elif [[ -z "${5}" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
elif [[ ! -z "${3}" ]] && [[ ! -f "${3}" ]]; then
	echo "Alternate database does not exist...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${5}" != "keep" ]] && [[ "${5}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 5th parameter...exiting"
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

"${shareScript}/clean_list.sh" "${1}"

# Set script directory
main_dir="${4}/GAMA_subs"
if [[ ! -d "${4}/GAMA_subs" ]]; then
	mkdir "${4}/GAMA_subs"
	mkdir "${4}/GAMA_subs/complete"
elif [[ ! -d  "${4}/GAMA_subs/complete" ]]; then
	mkdir "${4}/GAMA_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to check all isolates on the list against the newest ResGANNCBI DB
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" = "clobber" ]]; then
		rm ${processed}/${project}/${sample}/GAMA/${sample}_${ResGANNCBI_srst2_filename}.GAMA
	fi
	echo ${counter}
	# Check if counter is below max number of concurrent submissions
	if [ ${counter} -lt ${max_subs} ]; then
		# Check if the output file of GAMA exist, skip submission if so
		if [[ ! -f "${processed}/${project}/${sample}/GAMA/${sample}_${ResGANNCBI_srst2_filename}.GAMA" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -o GAMAAR_${sample}.out" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -e GAMAAR_${sample}.err" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -N GAMAAR_${sample}"   >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "cd ${shareScript}" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "\"${shareScript}/run_GAMA.sh\" \"${sample}\" \"${project}\" -c \"${3}\"" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_GAMAAR_complete.txt\"" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"

			#cd "${main_dir}"
			if [[ "${counter}" -lt "${last_index}" ]]; then
				qsub "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			else
				qsub -sync y "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			fi
			mv "${shareScript}/GAMAAR_${sample}.err" ${main_dir}
			mv "${shareScript}/GAMAAR_${sample}.out" ${main_dir}
		# Old data existed, skipping
		else
			echo -e $(date) > "${main_dir}/complete/${sample}_GAMAAR_complete.txt"
			echo "${project}/${sample} already has newest GAMA ResGANNCBI ${ResGANNCBI_srst2_filename}"
		fi
	# Counter is above max submission, must wait for previous ones to finish before moving on
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
		while :
		do
			# Check if timer is above max time allowed
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			# Check if waiting sample is finished
			if [ -f "${main_dir}/complete/${waiting_sample}_GAMAAR_complete.txt" ]; then
				# Check if current sample has etiher one of the output files from GAMA, skip analysis if so
				if [[ ! -f "${processed}/${project}/${sample}/GAMA/${sample}_${ResGANNCBI_srst2_filename}.GAMA" ]]; then
					echo  "Index is below max submissions, submitting"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -o GAMAAR_${sample}.out" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -e GAMAAR_${sample}.err" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -N GAMAAR_${sample}"   >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "cd ${shareScript}" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/run_GAMA.sh\" \"${sample}\" \"${project}\" -c \"${3}\"" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_GAMAAR_complete.txt\"" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"

					#cd "${main_dir}"
					if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					else
						qsub -sync y "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					fi
					mv "${shareScript}/GAMAAR_${sample}.err" ${main_dir}
					mv "${shareScript}/GAMAAR_${sample}.out" ${main_dir}
				# Old data existed, skipping
				else
					echo -e $(date) > "${main_dir}/complete/${sample}_GAMAAR_complete.txt"
					echo "${project}/${sample} already has newest GAMA ResGANNCBI ${ResGANNCBI_srst2_filename}"
				fi
				break
			# Wait 5 seconds and then check if "waiting" sample is complete
			else
				timer=$(( timer + 5 ))
				echo "sleeping for 5 seconds, so far slept for ${timer}"
				sleep 5
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_GAMAAR_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/GAMAAR_${waiting_sample}.out" ]]; then
			mv "${shareScript}/GAMAAR_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/GAMAAR_${waiting_sample}.err" ]]; then
			mv "${shareScript}/GAMAAR_${waiting_sample}.err" "${main_dir}"
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_GAMAAR_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/GAMAAR_${waiting_sample}.out" ]]; then
						mv "${shareScript}/GAMAAR_${waiting_sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/GAMAAR_${waiting_sample}.err" ]]; then
						mv "${shareScript}/GAMAAR_${waiting_sample}.err" "${main_dir}"
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

exit 0
