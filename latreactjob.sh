#!/bin/sh

#$ -M theavey@bu.edu
#$ -m beas
#$ -l h_rt=02:00:00
#$ -N latreact1

# name of gaussian input and output files (assumed to be in our home directory)
# EDIT THESE TWO LINES
# ARGS=\(0.5, 2.0, 2.0, 2., 2.0, 0.9, 0.75\) 100000 80 2000


# create a working area in local scratch space
SCRATCHDIR=/scratch/$USER
mkdir -p $SCRATCHDIR

# go there
cd $SCRATCHDIR

# copy input file into place
cp $HOME/lattice-reaction-model/*.py .

# set correct python environment from conda
# source activate py2.7
# this didn’t seem to work

# invoke gaussian
python runlat.py 1000 50 100 0.1 0.10 2.0 1.0 2.0 0.30 0.7


# print message
echo output is in /`hostname -s`$SCRATCHDIR