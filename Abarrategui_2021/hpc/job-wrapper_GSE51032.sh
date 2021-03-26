#!/bin/sh

# A simple job wrapper which logs the start date,
# loads the R module, then runs the R script 'job.R'.
# The end date is then logged.

echo "Start Date: "`date`
echo

module load R/4.0.3-foss-2020b
R -q -s -f GSE51032_FPR.R
R -q -s -f GSE51032_TPR.R

echo
echo "End Date:   "`date`
