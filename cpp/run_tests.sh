#!/bin/bash
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
BUILD_DIR="build"
TEST_TYPE="${TEST_TYPE:-all}"

echo -e "${BLUE}🧪 Running C++ Tests${NC}"
echo "========================"

# Check if build directory exists
if [ ! -d "$BUILD_DIR" ]; then
    echo -e "${YELLOW}⚠️  Build directory not found. Building first...${NC}"
    ./build_tests.sh
fi

cd "$BUILD_DIR"

# Set environment variables for optimal performance
export TBB_NUM_THREADS=$(nproc)
export OMP_NUM_THREADS=$(nproc)

# Run tests based on type
case "$TEST_TYPE" in
    "all")
        echo -e "${BLUE}🧪 Running all tests...${NC}"
        make run_tests
        ;;
    "multivector")
        echo -e "${BLUE}🧪 Running MultiVector tests...${NC}"
        ./multivector_tests
        ;;
    "distance")
        echo -e "${BLUE}🧪 Running DistanceMatrix tests...${NC}"
        ./distance_matrix_tests
        ;;
    "benchmarks")
        echo -e "${BLUE}⚡ Running benchmarks...${NC}"
        make run_benchmarks
        ;;
    "coverage")
        echo -e "${BLUE}📊 Running tests with coverage...${NC}"
        cmake .. -DCMAKE_BUILD_TYPE=Debug
        make -j$(nproc)
        ./multivector_tests
        ./distance_matrix_tests
        echo -e "${GREEN}✅ Coverage data generated in build directory${NC}"
        ;;
    *)
        echo -e "${RED}❌ Unknown test type: $TEST_TYPE${NC}"
        echo "Available types: all, multivector, distance, benchmarks, coverage"
        exit 1
        ;;
esac

echo -e "${GREEN}✅ Tests completed successfully!${NC}"

