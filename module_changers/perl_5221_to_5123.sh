#!/bin/sh -l

#$ -o regress_perl_5221.out
#$ -e regress_perl_5221.err
#$ -N regress_perl_5221
#$ -cwd
#$ -q short.q

# Unloads the python 3.6.1 module and replaces it with the 3.5.2 version, which all other modules function on
# Wouldnt unload normally so had to make this file...sorry
module unload perl/5.22.1
module load perl/5.12.3
