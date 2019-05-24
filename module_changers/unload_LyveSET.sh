#!/bin/sh -l

#$ -o unload_python_3.6.1.out
#$ -e unload_python_3.6.1.err
#$ -N unload_python_3.6.1
#$ -cwd
#$ -q short.q

# Script to load python 3.6.1 (only needed by busco currently)
# Wouldnt load normally so had to make this file...sorry
module unload Lyve-SET/1.1.4f
