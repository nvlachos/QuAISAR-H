#!/bin/sh -l

#$ -o amq-altc.out
#$ -e amq-altc.err
#$ -N amq-altc
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#List all currently loaded modules
#. ./module_changers/list_modules.sh

#
# Usage ./abl_mass_qsub_alt_csstar.sh path_to_list max_concurrent_submissions path_to_alt_database output_directory_for_scripts cloberness[keep|clobber] %ID(optional)[80|95|98|99|100]
#

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_alt_csstar.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions path_to_alt_database output_directory_for_scripts clobberness[keep|clobber] %ID(optional)[80|95|98|99|100]"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif ! [[ ${2} =~ $number ]] || [[ -z "${2}" ]]; then
	echo "${2} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${3}" ]]; then
	echo "alt_db parameter is empty...exiting"
	exit 1
elif [[ -z "${4}" ]]; then
	echo "No script output directory given...exiting"
	exit 4
elif [[ -z "${5}" ]] ; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${5}" != "keep" ]] && [[ "${5}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 5th parameter...exiting"
	exit 1
else
	clobberness="${5}"
fi

# Checks that value given for % Identity is one of the presets for csstar
if [[ "${6}" != 80 ]] && [[ "${6}" != 95 ]] && [[ "${6}" != 98 ]] && [[ "${6}" != 99 ]] && [[ "${6}" != 100 ]]; then
	echo "Identity is not one of the presets for csstar and therefore will fail, defaulting to 98..."
	sim="h"
	simnum=98
else
	if [ "${6}" == 98 ]; then
		sim="h"
	elif [ "${6}" == 80 ]; then
		sim="l"
	elif [ "${6}" == 99 ]; then
		sim="u"
	elif [ "${6}" == 95 ]; then
		sim="m"
	elif [ "${6}" == 100 ]; then
		sim="p"
	elif [ "${6}" == 40 ]; then
		sim="o"
	fi
	simnum=${6}
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [ "$line" ];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"


# Set counter and max submission variables
counter=0
max_subs=${2}

# Set script directory
main_dir="${4}/csstar_alt_subs"
if [[ ! -d "${4}/csstar_alt_subs" ]]; then
	mkdir "${4}/csstar_alt_subs"
	mkdir "${4}/csstar_alt_subs/complete"
elif [[ ! -d  "${4}/csstar_alt_subs/complete" ]]; then
	mkdir "${4}/csstar_alt_subs/complete"
fi

# format name being extracted from alt database
alt_database_path=$(basename -- "${3}")
echo "${alt_database_path}"
alt_database=$(echo ${alt_database_path##*/} | cut -d'.' -f1)
echo "${alt_database}"
alt_database=${alt_database//_srst2/}
echo "${alt_database}"


start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit qsub file for each isolate in the list, pausing if sunbissions reach max_concurrent_submissions until space frees up
while [ ${counter} -lt ${arr_size} ] ; do
	echo ${counter}
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	# Delete old results if clobbering
	if [[ "${5}" == "clobber" ]]; then
		rm ${processed}/${project}/${sample}/c-sstar/${sample}.${alt_database}.gapped_${simnum}_sstar_summary.txt
		rm -r ${processed}/${project}/${sample}/c-sstar/${alt_database}_gapped/
	fi

	# Check if there is an assembly to use
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		# Check if counter has reached max submissions yet
		if [[ ${counter} -lt ${max_subs} ]]; then
			echo  "Index is below max submissions, submitting"
			# Check if old results exist, and skip if so
			if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${alt_database}.gapped_${simnum}_sstar_summary.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				# Defaulting to gapped/98, change if you want to include user preferences
				echo -e "cd ${shareScript}" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_c-sstar_altDB.sh\" \"${sample}\" g "${sim}" \"${project}\" \"${3}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				#if [[ "${counter}" -lt "${last_index}" ]]; then
					qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
				#else
				#	if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
				#		qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
				#	else
				#		qsub -sync y "${main_dir}/csstn_${sample}_${start_time}.sh"
				#	fi
				#fi
			# Old data exists, skipping
			else
				echo "${project}/${sample} already has 0608"
				echo "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
			fi
		# Counter has exceeded max submissions limit
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			while :
			do
				# Check to see if max time limit has been reached, exit if so.
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check if the sample it is waiting on has been completed or if it has a usable assembly
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
					# Check if ol results exist and skip if so
					if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${alt_database}.gapped_${simnum}_sstar_summary.txt" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						# Defaulting to gapped/98, change if you want to include user preferences
						echo -e "cd ${shareScript}" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_c-sstar_altDB.sh\" \"${sample}\" g "${sim}" \"${project}\" \"${3}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						#if [[ "${counter}" -lt "${last_index}" ]]; then
							qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
						#else
						#	if [[ -d "${processed}/${project}/${sample}/c-sstar_plasmid" ]]; then
						#		qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
						#	else
						#		qsub -sync y "${main_dir}/csstn_${sample}_${start_time}.sh"
						#	fi
						#fi
					# Old data exists, skipping
					else
						echo "${project}/${sample} already has 0608"
						echo "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
					fi
					break
				# Wait 5 seconds and check if waiting sample is complete yet
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	# If no assembly file exists
	else
		echo "${project}/${sample} does not have an assembly file???"
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
ptimer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "${item} is complete normal"
		if [[ -f "${shareScript}/csstn_${waiting_sample}.out" ]]; then
			mv "${shareScript}/csstn_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/csstn_${waiting_sample}.err" ]]; then
			mv "${shareScript}/csstn_${waiting_sample}.err" "${main_dir}"
		fi
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 30 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]]; then
					echo "${item} is complete normal"
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
printf "%s %s" "abl_mass_qsub_alt_csstar.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_alt_csstar.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
