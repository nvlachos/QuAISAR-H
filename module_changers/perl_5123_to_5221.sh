#!/bin/sh -l

#$ -o advance_perl_5221.out
#$ -e advance_perl_5221.err
#$ -N advance_perl_5221
#$ -cwd
#$ -q short.q

# Unloads the perl 5.12.3 module and replaces it with the 5.22.1 version, which all other modules function on
# Wouldnt unload normally so had to make this file...sorry
module unload perl/5.12.3
module load perl/5.22.1
