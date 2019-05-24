#!/bin/sh -l

#$ -o list_mods.out
#$ -e list_mods.err
#$ -N list_mods
#$ -cwd
#$ -q short.q

# Script to print loaded modules to ensure the proper environment for running
# Wouldnt show mod list properly so made this file to run on cluster....sorry
module list
