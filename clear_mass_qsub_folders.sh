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
	rm -r blast16sID_*.out
	rm -r blast16sID_*.err
	rm -r blast16s_*.out
	rm -r blast16s_*.err
	rm -r ani_*.out
	rm -r ani_*.err
	rm -r ANI_*.out
	rm -r ANI_*.err
	rm -r BTQC_*.out
	rm -r BTQC_*.err
	rm -r BUSCO_*.out
	rm -r BUSCO_*.err
	rm -r getFASTQR1_*.out
	rm -r getFASTQR1_*.err
	rm -r getFASTQR2_*.out
	rm -r getFASTQR2_*.err
	rm -r csstn_*.out
	rm -r csstn_*.err
	rm -r csstp_*.out
	rm -r csstp_*.err
	rm -r kraka_*.out
	rm -r kraka_*.err
	rm -r krakr_*.out
	rm -r krakr_*.err
	rm -r gott_*.out
	rm -r gott_*.err
	rm -r mlst_*.out
	rm -r mlst_*.err
	rm -r pFinf_*.out
	rm -r pFinf_*.err
	rm -r pFinp_*.out
	rm -r pFinp_*.err
	rm -r SPAdn_*.out
	rm -r SPAdn_*.err
	rm -r SPAdp_*.out
	rm -r SPAdp_*.err
	rm -r plasFlow_*.out
	rm -r plasFlow_*.err
	rm -r pFinf_*.out
	rm -r pFinf_*.err
	rm -r PROKK_*.out
	rm -r PROKK_*.err
	rm -r PROKK_*.e*
	rm -r srst2AR_*.out
	rm -r srst2AR_*.err
	rm -r srst2MLST_*.out
	rm -r srst2MLST_*.err
	rm -r srst22MLST_*.out
	rm -r srst22MLST_*.err
	rm -r QUAST_*.out
	rm -r QUAST_*.err
	rm -r QC_*.out
	rm -r QC_*.err
	rm -r MLST_*.out
	rm -r MLST_*.err
	rm -r taxID_*.out
	rm -r taxID_*.err
	rm -r validate_*.out
	rm -r validate_*.err
	rm -r sum_*.out
	rm -r sum_*.err
	rm -r pFn_*.out
	rm -r pFn_*.err
	rm -r pFp_*.out
	rm -r pFp_*.err
fi
