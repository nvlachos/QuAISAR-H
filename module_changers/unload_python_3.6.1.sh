#!/bin/sh -l

#$ -o unload_python_3.6.1.out
#$ -e unload_python_3.6.1.err
#$ -N unload_python_3.6.1
#$ -cwd
#$ -q short.q

# Unloads the python 3.6.1 module and replaces it with the 3.5.2 version, which all other modules function on
# Wouldnt unload normally so had to make this file...sorry
module unload Python/3.6.1
module load Python/3.5.2
