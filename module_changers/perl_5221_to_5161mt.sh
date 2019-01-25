#!/bin/sh -l

#$ -o advance_perl_5221.out
#$ -e advance_perl_5221.err
#$ -N advance_perl_5221
#$ -cwd
#$ -q all.q

# Unloads the perl 5.22.1 module and replaces it with the 5.16.1-MT version, which all other modules function on
# Wouldnt unload normally so had to make this file...sorry
module unload perl/5.22.1
module load perl/5.16.1-MT
