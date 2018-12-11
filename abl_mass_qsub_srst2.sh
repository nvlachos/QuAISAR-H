#!/bin/sh -l

#$ -o ablmq_srst2.out
#$ -e ablmq_srst2.err
#$ -N ablmq_srst2
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./act_by_list.sh list_name(currently has to be placed in /scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR folder) description of list function
#
# script changes depending on what needs to be run through the list
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to act_by_list.sh, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./abl_mass_qsub_template.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_submissions"
	echo "Output location varies depending on which tasks are performed but will be found somewhere under ${processed}"
	exit 0
fi

arr=()
while IFS= read -r line || [[ "$line" ]];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

if [[ -z ${2} ]]; then
	max_subs=10
else
	max_subs=${2}
fi

# Loop through and act on each sample name in the passed/provided list
counter=0
main_dir="${share}/mass_subs/srst2_subs"
if [[ ! -d "${share}/mass_subs/srst2_subs" ]]; then
	mkdir "${share}/mass_subs/srst2_subs"
	mkdir "${share}/mass_subs/srst2_subs/complete"
elif [[ ! -d  "${share}/mass_subs/srst2_subs/complete" ]]; then
	mkdir "${share}/mass_subs/srst2_subs/complete"
fi

start_time=$(DATE)

while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	echo ${counter}
	if [ ${counter} -lt ${max_subs} ]; then
		if [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${resGANNOT_srst2_filename}_srst2__results.txt" ]] || [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${resGANNOT_srst2_filename}_srst2__results.txt" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module unload Python/2.7" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module unload Python/3.5.2" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module unload perl/5.22.1" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load Python/2.7.15" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load bowtie2/2.2.4" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load samtools/0.1.18" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load perl/5.16.1-MT" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "module load srst2" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
			echo -e "\"${shareScript}/run_srst2_on_singleDB.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			qsub -sync y "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			if [[ "${counter}" -lt "${last_index}" ]]; then
				qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			else
				qsub -sync y "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			fi
		else
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo "${project}/${sample} already has 20180608"
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
			if [ -f "${main_dir}/complete/${waiting_sample}_srst2_complete.txt" ]; then
				if [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__genes__${resGANNOT_srst2_filename}_srst2__results.txt" ]] && [[ ! -f "${processed}/${project}/${sample_name}/srst2/${sample_name}__fullgenes__${resGANNOT_srst2_filename}_srst2__results.txt" ]]; then
					echo "${waiting_sample} has completed, starting ${sample}"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module unload Python/2.7" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module unload Python/3.5.2" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module unload perl/5.22.1" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load Python/2.7.15" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load bowtie2/2.2.4" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load samtools/0.1.18" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load perl/5.16.1-MT" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "module load srst2" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
					echo -e "\"${shareScript}/run_srst2_on_singleDB.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					else
						qsub -sync y "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					fi
				else
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo "${project}/${sample} already has 20180608"
				fi
				break
			else
				timer=$(( timer + 5 ))
				echo "sleeping for 5 seconds, so far slept for ${timer}"
				sleep 5
			fi
		done
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
