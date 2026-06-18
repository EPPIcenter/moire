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
CMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE:-Release}"
PARALLEL_JOBS="${PARALLEL_JOBS:-$(nproc)}"

echo -e "${BLUE}🔨 Building C++ Tests${NC}"
echo "================================"

# Clean build directory if it exists
if [ -d "$BUILD_DIR" ]; then
    echo -e "${YELLOW}🧹 Cleaning build directory...${NC}"
    rm -rf "$BUILD_DIR"
fi

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Configure with CMake
echo -e "${BLUE}⚙️  Configuring with CMake...${NC}"
cmake .. \
    -DCMAKE_BUILD_TYPE="$CMAKE_BUILD_TYPE" \
    -DCMAKE_CXX_STANDARD=20 \
    -DCMAKE_CXX_STANDARD_REQUIRED=ON

# Build
echo -e "${BLUE}🔨 Building tests...${NC}"
make -j"$PARALLEL_JOBS"

echo -e "${GREEN}✅ Build completed successfully!${NC}"
echo ""
echo "Available executables:"
echo "  - multivector_tests"
echo "  - distance_matrix_tests"
echo "  - multivector_benchmarks"
echo "  - distance_matrix_benchmarks"
echo ""
echo "Run tests with:"
echo "  make run_tests"
echo "  make run_benchmarks"
echo ""
echo "Or run individual tests:"
echo "  ./multivector_tests"
echo "  ./distance_matrix_tests"

