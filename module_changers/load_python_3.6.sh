#!/bin/sh -l

#$ -o unload_python_3.6.out
#$ -e unload_python_3.6.err
#$ -N unload_python_3.6
#$ -cwd
#$ -q all.q

# Script to load python 3.6.1 (only needed by busco currently)
# Wouldnt load normally so had to make this file...sorry
module unload Python/3.5.2
module load Python/3.6.1
