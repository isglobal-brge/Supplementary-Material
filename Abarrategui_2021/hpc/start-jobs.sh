#!/bin/sh

#
# The error and log files are added to the 'logs' directory.
#

I_VALUES=`seq 1 100`
N_VALUES="20 30 40 50 60 70 80 90 100"

for I in ${I_VALUES}
do
    for N in ${N_VALUES}
    do
        
        sbatch --mem=15G --error=logs/'J'%j.err --output=logs/'J'%j.out --export=N=${N},I=${I} job-wrapper.sh
    done
done
    
