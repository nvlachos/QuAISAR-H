#!/bin/sh -l

#$ -o out_srst2.out
#$ -e out_srst2.err
#$ -N out_srst2
#$ -cwd
#$ -q all.q

# unloads the loaded modules of python, bowtie2, and samtools for srst2 that were necessary for srst2 to work correctly
# Wouldnt unload normally so had to make this file...sorry


module unload Python/2.7.11
module unload bowtie2/2.2.4
module unload samtools/0.1.18

module load Python/2.7
module load Python/3.5.2