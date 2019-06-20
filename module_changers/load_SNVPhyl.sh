#!/bin/sh -l

#$ -o load_snvphyl.out
#$ -e load_snvphy1.err
#$ -N load_snvphyl
#$ -cwd
#$ -q short.q

# Script to load python 3.6.1 (only needed by busco currently)
# Wouldnt load normally so had to make this file...sorry
#module unload Python/2.7
#module unload Python/3.5.2
module load snvphyl-galaxy-cli/1.3.0
module load Mash/2.0
