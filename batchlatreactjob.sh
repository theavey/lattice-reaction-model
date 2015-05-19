#!/bin/sh

# -pe omp 8

#$ -M theavey@bu.edu
#$ -m beas
#$ -l h_rt=02:00:00
#$ -N jobbatch

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
# this didn’t seem to work (neither above nor below)
# module load anaconda

# invoke gaussian
python runlat.py 10000 200 100 0.03 0.01 3.0 0.7 2.0 0.30 0.7 > out1 &
python runlat.py 10000 200 100 0.03 0.01 2.0 0.7 2.0 0.30 0.7 > out2 &
python runlat.py 10000 200 100 0.03 0.01 1.5 0.7 2.0 0.30 0.7 > out3 &
python runlat.py 10000 200 100 0.03 0.01 1.1 0.7 2.0 0.30 0.7 > out4 &
python runlat.py 10000 200 100 0.03 0.01 2.0 0.7 3.0 0.30 0.7 > out5 &
python runlat.py 10000 200 100 0.03 0.01 2.0 0.7 2.0 0.30 0.7 > out6 &
python runlat.py 10000 200 100 0.03 0.01 2.0 0.7 1.5 0.30 0.7 > out7 &
python runlat.py 10000 200 100 0.03 0.01 2.0 0.7 1.1 0.30 0.7 > out8 &
wait

# print message
echo output is in /`hostname -s`$SCRATCHDIR
