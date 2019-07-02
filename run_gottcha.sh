#!/bin/sh -l

#$ -o run_gottcha.out
#$ -e run_gottcha.err
#$ -N run_gottcha
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Runs the gottcha classification tool (now only at species level) which identifies the most likely taxonomic classification for the sample
#
# Usage ./run_gottcha.sh sample_name run_id
#
# requires modules gottcha and perl 5.12.3
#
# !Version 1
#

ml perl/5.12.3  gottcha

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_gottcha.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_gottcha.sh   sample_name    run_id"
	echo "Output is saved to ${processed}/miseq_run_id/sample_name/gottcha/"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id name supplied to run_gottcha.sh, exiting"
	exit 1
fi

# Sets the output folder of gottcha classifier to the gottcha folder under the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}"

# Create necessary output directories
echo "Running gottcha Taxonomic Classifier"
if [ ! -d "$OUTDATADIR/gottcha" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/gottcha"
	mkdir -p "$OUTDATADIR/gottcha"
fi
if [ ! -d "$OUTDATADIR/gottcha/gottcha_S" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/gottcha/gottcha_S"
	mkdir -p "$OUTDATADIR/gottcha/gottcha_S"
fi

# Needs perl v5.12.3 to function properly
. "${shareScript}/module_changers/perl_5221_to_5123.sh"

### Gottcha Taxonomy Classifier ### in species mode
gottcha.pl --mode all --outdir "${OUTDATADIR}/gottcha/gottcha_S" --input "${OUTDATADIR}/trimmed/${1}.paired.fq" --database "${gottcha_db}"

# Return perl to 5.22.1
. "${shareScript}/module_changers/perl_5123_to_5221.sh"

# Create the krona graphs from each of the analyses
ktImportText "${OUTDATADIR}/gottcha/gottcha_S/${1}_temp/${1}.lineage.tsv" -o "${OUTDATADIR}/gottcha/${1}_species.krona.html"

#Create a best hit from gottcha1 file
"${shareScript}/best_hit_from_gottcha1.sh" "${1}" "${2}"

ml -perl/5.12.3 -gottcha

#Script exited gracefully (unless something else inside failed)
exit 0
