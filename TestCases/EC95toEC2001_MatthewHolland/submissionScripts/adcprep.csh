#! /bin/csh
#BSUB -J adcprep
#BSUB -o adcprep.%J
#BSUB -e adcprep.%J
#BSUB -W 15
#BSUB -n 1
#BSUB -q ccee

./adcprep --np 256 --partmesh
./adcprep --np 256 --prepall
