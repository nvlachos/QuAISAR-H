#!/bin/sh -l

#$ -o run_ANI.out
#$ -e run_ANI.err
#$ -N run_ANI
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh
#   ${mod_changers}/list_modules.sh

ml Python3/3.5.2 pyani/0.2.7 mashtree/0.29

#
# Script to calculate the average nucleotide identity of a sample to numerous other samples from the same genus (genus dependent)
# The most similar match is identified and provided for confirmation
#
# Usage ./run_ANI.sh sample_name   DB(for looking up reference, just relative path, also is genus)   Species   run_id  list_samples_to_include(optional)
#
# Python/3.5.2 (pyani is located in Nick_DIR/script folder, not run from scicomp module)
#



# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_ANI.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_ANI.sh sample_name ani_database(which is also genus) species run_id list_of_samples_to_include(optional)"
	echo "Output is saved to in ${processed}/sample_name/ANI"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty database name supplied to run_ANI.sh. Second argument should be a genus found in ${local_DBs}/ANI/  ...Exiting"
	exit 1
elif [ ! -s "${local_DBs}/aniDB/${2,}" ]; then
	echo "The genus does not exist in the ANI database. This will be noted and the curator of the database will be notified. However, since nothing can be done at the moment....exiting"
	# Create a dummy folder to put non-results into (if it doesnt exist
	if [ ! -d "${processed}/${4}/${1}/ANI" ]; then  #create outdir if absent
		echo "${processed}/${4}/${1}/ANI"
		mkdir -p "${processed}/${4}/${1}/ANI"
	fi
	# Write non-results to a file in ANI folder
	echo "No matching ANI database found for ${2}(genus)" >> "${processed}/${4}/${1}/ANI/best_ANI_hits_ordered(${1}_vs_${2}).txt"
	# Add genus to list to download and to database
	global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
	echo "ANI: ${2} - Found as ${1} on ${global_time}" >> "${shareScript}/maintenance_To_Do.txt"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty species name supplied to run_ANI.sh. Third argument should be the suspected species of the sample. Exiting"
	exit 1
elif [ -z "$4" ]; then
	echo "Empty miseq_run_id name supplied to run_ANI.sh. Fourth argument should be the run id. Exiting"
	exit 1
elif [ ! -z "$5" ]; then
	others="true"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Started ANI at ${start_time}"

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${4}/${1}"

# Checks to see if an ANI folder already exists and creates it if not
if [ ! -d "$OUTDATADIR/ANI" ]; then
	echo "No $OUTDATADIR/ANI directory, exiting"
	exit
fi

# Gets persons name to use as email during entrez request to identify best matching sample
me=$(whoami)
#echo ${me}"___"${1}___${2}___${3}___${4}

# Sets the genus as the database that was passed in (The $2 seemed to be getting changed somewhere, so I just set it as a local variable)
genus_in=${2}

#Extracts all %id against the query sample (excluding itself) and writes them to file
if [[ ! -d "${OUTDATADIR}/ANI/localANIDB" ]]; then
	mkdir "${OUTDATADIR}/ANI/localANIDB"
	for (( i=0; i<n; i++ ));
	do
		temp_ref=$(find ${local_DBs}/aniDB/${genus_in,,} -maxdepth 1 -type f -name "*${samples[i]}.fna.gz")
		echo "Trying to copy ${temp_ref} --- ${samples[i]}.fna.gz"
		if [[ -f ${temp_ref} ]]; then
			cp "${temp_ref}" "${OUTDATADIR}/ANI/localANIDB"
		else
			echo "Could not copy ${temp_ref} (*${samples[i]}.fna.gz)"
		fi
	done
	gunzip "${OUTDATADIR}/ANI/localANIDB/"*
else
	echo "Already/still has its localANIDB folder"
fi

rename 's/.fna$/.fasta/' ${OUTDATADIR}/ANI/localANIDB/*.fna

mashtree --numcpus ${procs} *.fna --tempdir ${OUTDATADIR}/ANI/temp > ${OUTDATADIR}/ANI/"${genus_in}_and_${1}_mashtree.dnd";

rename 's/.fasta$/.fna/' ${OUTDATADIR}/ANI/localANIDB/*.fna

end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "ENDed ANI at ${end_time}"

#Script exited gracefully (unless something else inside failed)
exit 0
