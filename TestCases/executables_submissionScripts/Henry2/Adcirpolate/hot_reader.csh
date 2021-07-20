#! /bin/csh
#BSUB -J adcirpolate
#BSUB -o adcirpolate.%J
#BSUB -e adcirpolate.%J
#BSUB -W 2880
#BSUB -n 256
#BSUB -q ccee

mpirun ./adcirpolate
