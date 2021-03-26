Project Demo
------------

- README.txt

  This text!

- inputs

  Directory containing input files used by jobs.

- outputs

  Directory containing output files used by jobs.

- logs

  Diectory containing error and log files generated for jobs.

- start-jobs.sh

  A shell script which start a set of jobs.

- job-wrapper.sh

  A shell script which runs the R script which perform the job.

- job.R

  A R script which contains the job to be executed.

Notes
-----

- squeue -u username

  The command "squeue -u username" will list the jobs in the queue
  for the specified user.

