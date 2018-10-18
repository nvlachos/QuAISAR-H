#!/bin/sh -l

#$ -o unload_snvphyl.out
#$ -e unload_snvphy1.err
#$ -N unload_snvphyl
#$ -cwd
#$ -q all.q

# Script to load python 3.6.1 (only needed by busco currently)
# Wouldnt load normally so had to make this file...sorry
module unload snvphyl-galaxy-cli/1.3.0
module load Python/2.7
module load Python/3.5.2

