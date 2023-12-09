#!/bin/bash
  
#SBATCH --account=nlscience
#SBATCH --chdir=./ # Working directory
#SBATCH --mail-user=battistn@tcnj.edu # Who to send emails to
#SBATCH --mail-type=ALL             # Send emails on start, end and failure
#SBATCH --job-name=IBFE_5b          # Job name
#SBATCH --output=job.%j.out         # Name of stdout output file (%j expands to jobId)
#SBATCH --nodes=1                   # Total number of nodes (a.k.a. servers) requested
#SBATCH --ntasks=1                  # Total number of mpi tasks requested
#SBATCH --partition=normal          # Partition (a.k.a.queue) to use
#SBATCH --requeue                   # restart flag if HPC goes down

###SBATCH --partition=nolimit         # Partition (a.k.a.queue) to use
###SBATCH --time=60-00:00:00          # Run time (days-hh:mm:ss)

# Launch a serial job
echo "Starting @ "`date`
./main2d input2d
echo "Completed @ "`date`
