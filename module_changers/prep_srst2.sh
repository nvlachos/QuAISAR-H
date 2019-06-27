#!/bin/sh -l

#$ -o out_srst2.out
#$ -e out_srst2.err
#$ -N out_srst2
#$ -cwd
#$ -q short.q

# unloads the loaded modules of python, bowtie2, and samtools for srst2 that were necessary for srst2 to work correctly
# Wouldnt unload normally so had to make this file...sorry

#module unload Python/2.7
#module unload python/3.5.2
#module unload perl/5.22.1

#module load Python/2.7.15
#module load bowtie2/2.2.4
#module load samtools/0.1.18
#module load perl/5.16.1-MT
#module load srst2


ml -Python/2.7 -Python/3.5.2 -perl/5.22.1
ml Python/2.7.11 bowtie2/2.2.4 samtools/0.1.18 perl/5.16.1-MT srst2
