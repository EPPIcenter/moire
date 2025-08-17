#!/bin/bash
set -e

# Set TBB thread count to match the number of physical cores
# This helps avoid oversubscription which can hurt performance
export TBB_NUM_THREADS=$(nproc)

# Compile the benchmark test with TBB support and C++20
# Use -O2 instead of -O3 for more realistic benchmarking
# Add -fno-inline to prevent aggressive inlining
# g++ -std=c++20 -O2 -fno-inline -I. -I../src -o multivector_benchmark multivector_benchmark.cpp -pthread -ltbb
g++ -std=c++20 -O3 -march=native -I. -I../src -o multivector_benchmark multivector_benchmark.cpp -pthread -ltbb

# Run the benchmark
#./multivector_benchmark 
