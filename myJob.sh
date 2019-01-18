#!/bin/bash
#SBATCH -J myMPI # Job Name
#SBATCH -o myMPI.o%j # Name of the output file
#SBATCH -N 8   # number of nodes requested 
#SBATCH -n 16  # total number of mpi tasks requested
#SBATCH -p normal  # queue (partition) -- normal, development, etc.
#SBATCH -t 1:00:00  # Run time (hh:mm:ss) - 1 hours

set -x # Echo commands 

ibrun ./Helloworld # Run the MPI executable
