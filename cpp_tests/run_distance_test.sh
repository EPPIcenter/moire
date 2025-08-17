#!/bin/bash
set -e

# Set TBB thread count to match the number of physical cores
# This helps avoid oversubscription which can hurt performance
export TBB_NUM_THREADS=$(nproc)

g++ -std=c++20 -O3 -march=native -I. -I../src -o distance_test test_distance_matrix.cpp -pthread -ltbb

# Run the benchmark
./distance_test 
