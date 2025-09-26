#!/bin/bash -x
# submit a chain of jobs with dependency
# number of jobs to submit
NO_OF_JOBS=4
# define jobscript
JOB_SCRIPT=runDALES.sh
echo "sbatch ${JOB_SCRIPT}"
JOBID=$(sbatch ${JOB_SCRIPT} 2>&1 | awk '{print $(NF)}')
I=2
while [ ${I} -le ${NO_OF_JOBS} ]; do
echo "sbatch --dependency=afterok:${JOBID} ${JOB_SCRIPT}"
JOBID=$(sbatch --dependency=afterok:${JOBID} ${JOB_SCRIPT} 2>&1 | awk '{print $(NF)}')
let I=${I}+1
done
