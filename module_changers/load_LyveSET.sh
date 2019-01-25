#!/bin/sh -l

#$ -o load_python_3.6.1.out
#$ -e load_python_3.6.1.err
#$ -N load_python_3.6.1
#$ -cwd
#$ -q all.q

# Script to load python 3.6.1 (only needed by busco currently)
# Wouldnt load normally so had to make this file...sorry
module load Lyve-SET/1.1.4f
