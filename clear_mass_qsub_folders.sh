#!/bin/sh -l

#$ -o clear_mqs.out
#$ -e clear_mqs.err
#$ -N clear_mqs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
. "${mod_changers}/pipeline_mods"

#
# Usage ./clear_mass_qsub_fodlers [1,2,3] (1-qsub folders, 2-qsub outs/errs, 3-Both)
#
# script changes naming structure of SPAdes output to include isolate name for every contig and removes coverage info
#
if [[ ${1} -eq 1 ]] || [[ ${1} -eq 3 ]]; then
	for folder in ${share}/mass_subs/*
	do
		if [[ -d ${folder} ]]; then
			for folder1 in ${folder}
			do
				echo "Deleting every sh script in ${folder1}"
				rm -r ${folder1}/*.sh
				echo "Deleting every sh script in ${folder1}/complete"
				rm -r ${folder1}/complete/*.txt
			done
		fi
	done
fi

if [[ ${1} -eq 2 ]] || [[ ${1} -eq 3 ]]; then
	rm  blast16sID_*.out
	rm  blast16sID_*.err
	rm  blast16s_*.out
	rm  blast16s_*.err
	rm  ani_*.out
	rm  ani_*.err
	rm  ANI_*.out
	rm  ANI_*.err
	rm  BTQC_*.out
	rm  BTQC_*.err
	rm  BUSCO_*.out
	rm  BUSCO_*.err
	rm  getFASTQR1_*.out
	rm  getFASTQR1_*.err
	rm  getFASTQR2_*.out
	rm  getFASTQR2_*.err
	rm  csstn_*.out
	rm  csstn_*.err
	rm  csstp_*.out
	rm  csstp_*.err
	rm  kraka_*.out
	rm  kraka_*.err
	rm  krakr_*.out
	rm  krakr_*.err
	rm  gott_*.out
	rm  gott_*.err
	rm  mlst_*.out
	rm  mlst_*.err
	rm  pFinf_*.out
	rm  pFinf_*.err
	rm  pFinp_*.out
	rm  pFinp_*.err
	rm  SPAdn_*.out
	rm  SPAdn_*.err
	rm  SPAdp_*.out
	rm  SPAdp_*.err
	rm  plasFlow_*.out
	rm  plasFlow_*.err
	rm  pFinf_*.out
	rm  pFinf_*.err
	rm  PROKK_*.out
	rm  PROKK_*.err
	rm  PROKK_*.e*
	rm  srst2AR_*.out
	rm  srst2AR_*.err
	rm  srst2MLST_*.out
	rm  srst2MLST_*.err
	rm  srst22MLST_*.out
	rm  srst22MLST_*.err
	rm  QUAST_*.out
	rm  QUAST_*.err
	rm  QC_*.out
	rm  QC_*.err
	rm  MLST_*.out
	rm  MLST_*.err
	rm  taxID_*.out
	rm  taxID_*.err
	rm  validate_*.out
	rm  validate_*.err
	rm  sum_*.out
	rm  sum_*.err
	rm  pFn_*.out
	rm  pFn_*.err
	rm  pFp_*.out
	rm  pFp_*.err
	rm aniB_*.out
	rm aniB_*.err
	rm aniM_*.out
	rm aniM_*.err
	rm core.*
fi
