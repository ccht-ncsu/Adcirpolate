#! /bin/csh
#BSUB -J adcprepHot
#BSUB -o adcprep.%J
#BSUB -e adcprep.%J
#BSUB -W 15
#BSUB -n 1
#BSUB -q ccee

./adcprep <in.prephot
