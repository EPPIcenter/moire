#!/bin/bash
set -e

# Set TBB thread count to match the number of physical cores
# This helps avoid oversubscription which can hurt performance
export TBB_NUM_THREADS=$(nproc)

# Compile the benchmark test with TBB support and C++20
# Use -O2 instead of -O3 for more realistic benchmarking
# Add -fno-inline to prevent aggressive inlining
# g++ -std=c++20 -O2 -fno-inline -I/home/mmurphy/R/x86_64-pc-linux-gnu-library/4.4/RcppParallel/include/tbb -I. -I../src -L/home/mmurphy/R/x86_64-pc-linux-gnu-library/4.4/RcppParallel/lib -o multivector_test multivector_test.cpp -ltbb
# g++ -std=c++20 -O2 -fno-inline -I/home/mmurphy/R/x86_64-pc-linux-gnu-library/4.4/RcppParallel/include -L/home/mmurphy/R/x86_64-pc-linux-gnu-library/4.4/RcppParallel/lib -I. -I../src -o multivector_test multivector_test.cpp -ltbb
g++ -std=c++20 -O2 -fno-inline -fopenmp -I. -I../src -o multivector_test multivector_test.cpp -ltbb
#clang++ -std=c++20 -O2 -fno-inline -fopenmp -I. -I../src -o multivector_test multivector_test.cpp -ltbb

# Run the test
./multivector_test
