#!/bin/sh

#$ -pe omp 8

#$ -M theavey@bu.edu
#$ -m beas
#$ -l h_rt=12:00:00
#$ -N tempdepend6


# create a working area in local scratch space
SCRATCHDIR=/scratch/$USER
mkdir -p $SCRATCHDIR

# go there
cd $SCRATCHDIR

# copy code files into place
cp $HOME/lattice-reaction-model/*.py .

# invoke python scripts to be run in parallel
python runlat.py 35000 120 500 0.02 0.375 2.0 0.5 0.5 0.30 0.1 > out1 &
python runlat.py 35000 120 500 0.02 0.375 2.0 0.5 0.5 0.30 0.1 > out2 &
python runlat.py 35000 120 500 0.02 0.375 2.0 0.5 0.5 0.30 0.5 > out3 &
python runlat.py 35000 120 500 0.02 0.375 2.0 0.5 0.5 0.30 0.5 > out4 &
python runlat.py 35000 120 500 0.02 0.375 2.0 0.5 0.5 0.30 1.0 > out5 &
python runlat.py 35000 120 500 0.02 0.375 2.0 0.5 0.5 0.30 1.0 > out6 &
python runlat.py 35000 120 500 0.02 0.375 2.0 0.5 0.5 0.30 3.0 > out7 &
python runlat.py 35000 120 500 0.02 0.375 2.0 0.5 0.5 0.30 3.0 > out8 &
wait

# print message
echo output is in /net/scc-`hostname -s`$SCRATCHDIR
