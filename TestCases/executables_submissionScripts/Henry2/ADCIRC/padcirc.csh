#! /bin/csh
#BSUB -J padcirc
#BSUB -o padcirc.%J
#BSUB -e padcirc.%J
#BSUB -W 2880
#BSUB -n 256
#BSUB -q ccee

mpirun ./padcirc
