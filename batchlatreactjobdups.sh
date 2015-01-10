#!/bin/sh

# -pe omp 8

#$ -M theavey@bu.edu
#$ -m beas
#$ -l h_rt=04:00:00
#$ -N jobbatch4

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
python runlat.py 10000 200 100 0.03 0.01 4.0 0.7 0.25 0.30 0.1 > out1 &
python runlat.py 10000 200 100 0.03 0.01 4.0 0.7 0.25 0.30 0.1 > out2 &
python runlat.py 10000 200 100 0.03 0.01 2.0 0.7 0.50 0.30 0.1 > out3 &
python runlat.py 10000 200 100 0.03 0.01 2.0 0.7 0.50 0.30 0.1 > out4 &
python runlat.py 10000 200 100 0.03 0.01 0.50 0.7 2.0 0.30 0.1 > out5 &
python runlat.py 10000 200 100 0.03 0.01 0.50 0.7 2.0 0.30 0.1 > out6 &
python runlat.py 10000 200 100 0.03 0.01 0.25 0.7 4.0 0.30 0.1 > out7 &
python runlat.py 10000 200 100 0.03 0.01 0.25 0.7 4.0 0.30 0.1 > out8 &
wait

# print message
echo output is in /`hostname -s`$SCRATCHDIR
