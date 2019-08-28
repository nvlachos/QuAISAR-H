#!/bin/sh -l

#$ -o run_GAMA.out
#$ -e run_GAMA.err
#$ -N run_GAMA
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Runs the GAMA AR classification tool
#
# Usage ./run_GAMA.sh sample_name run_id -c|p [path_to_alt_DB]
#
# requires modules blat Python/2.7.3
#
# !Version 1
#

ml blat Python/2.7.3

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_GAMA.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_GAMA.sh   sample_name    run_id	c|p	[path_to_alt_DB]"
	echo "Output is saved to ${processed}/miseq_run_id/sample_name/GAMA/"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id name supplied to run_GAMA.sh, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty assembly source supplied to run_GAMA.sh (must be -c or -p, chromosome or plasmid respectively), exiting"
	exit 1
elif [[ "${3}" != "-c" &&  "${3}" != "-p" ]]; then
	echo "Incorrect assembly source supplied to run_GAMA.sh (must be -c or -p, chromosome or plasmid respectively), exiting"
	exit 1
elif [ ! -z "${4}" ]; then
	ARDB="${4}"
else
	ARDB="${ResGANNCBI_srst2}"
fi

# Sets the output folder of GAMA classifier to the GAMA folder under the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}"

# Create necessary output directories
echo "Running GAMA Taxonomic Classifier"
if [ ! -d "$OUTDATADIR/GAMA" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/GAMA"
	mkdir -p "$OUTDATADIR/GAMA"
fi

OUTDATA="${OUTDATADIR}"

if [[ "${3}" == "-c" ]]; then
	assembly_source="${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta"
	OUTDATADIR="${OUTDATADIR}/GAMA"
elif [[ "${3}" == "-p" ]]; then
	assembly_source="${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta"
	OUTDATADIR="${OUTDATADIR}/GAMA_plasFlow"
else
	echo "Unknown Assembly source identifier, exiting"
	exit 5564
fi
### GAMA AR Classifier ### in species mode
python2 GAMA_4.6_ResGANNOT_SciComp_Exe.py "${OUTDATA}/Assembly/${1}_scaffolds_trimmed.fasta" "${ARDB}" "${OUTDATADIR}/${1}_${ResGANNOT_srst2_filename}.GAMA"

ml -blat -Python/2.7.3

#Script exited gracefully (unless something else inside failed)
exit 0
