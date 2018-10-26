#!/bin/sh -l

#$ -o run_coreestimate.out
#$ -e run_coreestimate.err
#$ -N run_coreestimate
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import list of mods used during pipeline analysis (or downstream)
. "${mod_changers}/pipeline_mods"

#
# Runs Toms core genome estimator from BAM file on a Lyve-SET run
#
# Usage ./run_coreestimate.sh group_name (must have sample list file in Phylogeny_analyses folder before starting as group_name.samples. First sample on list will be reference)
#
# No modules required yet (tested with perl 5.12.3 being used)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_coreestimate.sh, exiting"
	exit 1
elif [[ -z "${1}" ]] || [ ! -s "${share}/Phylogeny_analyses/${1}/${1}.samples" ]; then
	echo "Empty group name supplied to run_core_estimate.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_coreestimate.sh group_name (to get core genome size of all samples)"
	echo "Phylogeny_analyses folder must contain list of samples to be compared labeled by group_name. First sample on list will be reference"
	echo "Output is saved to ${share}/Phylogeny_analyses/group_name/Lyve-SET"
	exit 0
fi

ref=$(head -n 1 "${share}/Phylogeny_analyses/${1}/${1}.samples")
ref=$(echo "${ref}" | cut -d'/' -f2 | tr -d '[:space:]')
#owd=$(pwd)
if [[ "${2}" = "-LyveSET" ]] || [[ "${2}" = "-L" ]] || [[ "${2}" = "-l" ]]; then
	if [[ -s "${share}/Phylogeny_analyses/${1}/LyveSET/reference/reference.fasta" ]]; then 
		#cd "${share}/Phylogeny_analyses/${1}/LyveSET/output"
		perl "${shareScript}/estimate_core_genome_from_bam.pl" -bam "${share}/Phylogeny_analyses/${1}/LyveSET/bam" -genome "${share}/Phylogeny_analyses/${1}/LyveSET/reference/reference.fasta" -out "${share}/Phylogeny_analyses/${1}/LyveSET/output/" -depth 10 2 #> "${share}/Phylogeny_analyses/${1}/LyveSET/output/genome_core.txt"
	else
		echo "Reference FASTA does not exist...exiting"
	fi
elif [[ "${2}" = "-SNVPhyl" ]] || [[ "${2}" = "-S" ]]  || [[ "${2}" = "-s" ]]; then
	reference=$(find ${share}/Phylogeny_analyses/${1}/SNVPhyl -name "reference*.fasta" -type f)
	if [[ ! -z  "${reference}" ]]; then
		#cd "${share}/Phylogeny_analyses/${1}/SNVPhyl/output"
		perl "${shareScript}/estimate_core_genome_from_bam.pl" -bam "${share}/Phylogeny_analyses/${1}/SNVPhyl/output/bam" -genome "${share}/Phylogeny_analyses/${1}/SNVPhyl/${reference}" -out "${share}/Phylogeny_analyses/${1}/SNVPhyl/output/" -depth 10 2 #> "${share}/Phylogeny_analyses/${1}/SNVPhyl/output/genome_core.txt"
	else
		echo "Reference FASTA does not exist...exiting"
	fi
fi 
cd ${owd}
#Script exited gracefully (unless something else inside failed)
exit 0
