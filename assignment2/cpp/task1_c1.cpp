#!/bin/bash
#SBATCH -J myMPI           # job name
#SBATCH -o myMPI.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 68              # total number of mpi tasks requested
#SBATCH -p normal     # queue (partition) -- normal, development, etc.
#SBATCH -t 07:30:00        # run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=hmayani@cs.stonybrook.edu
#SBATCH --mail-type=begin  # email me when the job starts
#SBATCH --mail-type=end    # email me when the job finishes

icpc -std=c++11 ParRadixCount.cpp -o radixC
./radixC 13 1
./radixC 14 1
./radixC 15 1
./radixC 16 1
./radixC 17 1
./radixC 18 1
./radixC 19 1
./radixC 20 1

icpc -std=c++11 ParQuickSort.cpp -o quick
./quick 13 13 1 
./quick 14 13 1
./quick 15 13 1
./quick 16 13 1
./quick 17 13 1
./quick 18 13 1
./quick 19 13 1
./quick 20 13 1

