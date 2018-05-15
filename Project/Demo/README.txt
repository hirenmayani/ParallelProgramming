module load intel/2016.3.210 
icpc cilk_histo.cpp -std=c++11 -o cilk_hist.out
icpc cilk_wc.cpp -std=c++11 -o cilk_wc.out
icpc -qopenmp openmp_wc.cpp -std=c++11 -o omp_wc.out
icpc -qopenmp openmp_hist.c -std=c++11 -o omp_hist.out
