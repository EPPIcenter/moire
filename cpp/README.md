# C++ Testing Framework for Moire

This directory contains a comprehensive C++ testing framework for the Moire package, designed to provide robust testing, benchmarking, and validation of C++ components.

## 🏗️ **Architecture**

### **Test Framework Components**

1. **Custom Test Framework** (`test_framework.hpp`)
   - Lightweight, header-only testing framework
   - Comprehensive assertion macros
   - Test result tracking and reporting
   - Performance monitoring

2. **CMake Build System** (`CMakeLists.txt`)
   - Modern CMake configuration
   - Automatic dependency management
   - Cross-platform support
   - Coverage and profiling support

3. **Test Suites**
   - `MultiVectorTestSuite`: Tests for MultiVector class
   - `DistanceMatrixTestSuite`: Tests for DistanceMatrix class
   - Extensible framework for additional test suites

4. **Build Scripts**
   - `build_tests.sh`: Automated build process
   - `run_tests.sh`: Test execution with various options

## 🚀 **Quick Start**

### **Prerequisites**
- C++20 compatible compiler (GCC 10+, Clang 12+, MSVC 2019+)
- CMake 3.20+
- TBB (Threading Building Blocks)
- OpenMP
- Catch2 (automatically downloaded)
- Google Benchmark (automatically downloaded)

### **Building Tests**
```bash
# Build all tests
./build_tests.sh

# Or use CMake directly
mkdir build && cd build
cmake ..
make -j$(nproc)
```

### **Running Tests**
```bash
# Run all tests
./run_tests.sh

# Run specific test types
TEST_TYPE=multivector ./run_tests.sh
TEST_TYPE=distance ./run_tests.sh
TEST_TYPE=benchmarks ./run_tests.sh
TEST_TYPE=coverage ./run_tests.sh
```

### **Using Justfile Integration**
```bash
# Build C++ tests
just cpp-build

# Run C++ tests
just cpp-test

# Run benchmarks
just cpp-benchmark

# Run with coverage
just cpp-coverage

# Run specific test suite
just cpp-test-suite multivector
```

## 🧪 **Test Framework Features**

### **Assertion Macros**
```cpp
ASSERT_TRUE(condition)           // Assert condition is true
ASSERT_FALSE(condition)          // Assert condition is false
ASSERT_EQ(expected, actual)     // Assert equality
ASSERT_NE(expected, actual)     // Assert inequality
ASSERT_LT(left, right)          // Assert less than
ASSERT_LE(left, right)          // Assert less than or equal
ASSERT_GT(left, right)          // Assert greater than
ASSERT_GE(left, right)          // Assert greater than or equal
ASSERT_THROWS(statement)         // Assert statement throws
ASSERT_NO_THROW(statement)      // Assert statement doesn't throw
```

### **Test Suite Structure**
```cpp
class MyTestSuite : public TestSuite {
public:
    MyTestSuite() : TestSuite("MyTestSuite") {
        add_test("Test Name", [this]() { test_function(); });
        add_test("Another Test", [this]() { another_test(); });
    }

private:
    void test_function() {
        // Test implementation
        ASSERT_TRUE(some_condition);
    }
};
```

### **Performance Testing**
```cpp
void test_performance() {
    auto start = std::chrono::high_resolution_clock::now();
    // Code to test
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    
    ASSERT_LT(duration, 1000); // Should complete in less than 1 second
}
```

## 📊 **Test Categories**

### **1. Unit Tests**
- **Construction Tests**: Verify object creation and initialization
- **Access Tests**: Test element access and modification
- **Copy/Move Tests**: Test copy constructors, assignment operators
- **Edge Cases**: Test boundary conditions and error cases

### **2. Integration Tests**
- **Component Interaction**: Test how different components work together
- **Data Flow**: Test data passing between components
- **Error Propagation**: Test error handling across components

### **3. Performance Tests**
- **Benchmarking**: Measure execution time and memory usage
- **Scalability**: Test performance with different data sizes
- **Optimization**: Verify performance improvements

### **4. Validation Tests**
- **Mathematical Properties**: Verify distance matrix properties
- **Numerical Stability**: Test with edge case values
- **Precision**: Test floating-point accuracy

## 🔧 **Configuration**

### **Environment Variables**
```bash
export CMAKE_BUILD_TYPE=Release    # Build type (Debug/Release)
export PARALLEL_JOBS=8            # Number of parallel build jobs
export TBB_NUM_THREADS=4          # TBB thread count
export OMP_NUM_THREADS=4          # OpenMP thread count
export TEST_TYPE=all              # Test type to run
```

### **CMake Options**
```bash
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_CXX_STANDARD=20 \
    -DENABLE_COVERAGE=ON \
    -DENABLE_PROFILING=ON
```

## 📈 **Benchmarking**

### **Running Benchmarks**
```bash
# Run all benchmarks
TEST_TYPE=benchmarks ./run_tests.sh

# Run specific benchmark
./build/multivector_benchmarks1
./build/multivector_benchmarks2 100 100 100
./build/multivector_regression_benchmarks

# Emit machine-readable CSV for CI / regression tracking
MOIRE_BENCHMARK_CSV=1 ./build/multivector_regression_benchmarks

# prob_any_missing (inclusion-exclusion) sweep: total_alleles x coi, vectorized + operator (Gray-code baseline)
# Run validates Gray vs combination (same formula) within tol 1e-5 before timing.
./build/prob_any_missing_benchmarks
MOIRE_BENCHMARK_CSV=1 ./build/prob_any_missing_benchmarks  # regression CSV; reps/iter adaptive by size

# Custom prob_any_missing: pass event probs and numEvents
./build/prob_any_missing_benchmarks 0.1 0.2 0.3 0.4 10
```

### **Benchmark Features**
- **Multiple Iterations**: Configurable warmup and measurement runs
- **Regression suite**: `multivector_regression_benchmarks` uses fixed seed (42) for repeatable size/operation sweeps (reduce, sum, logsumexp, transform; seq vs parallel).
- **CSV output**: Set `MOIRE_BENCHMARK_CSV=1` to print CSV for baseline capture and CI trend checks. Multivector: `scenario,dim0,dim1,dim2,operation,policy,time_ms`. prob_any_missing: `scenario,total_alleles,coi,operation,time_ms` (time_ms is per-call).
- **P(any missing) accuracy**: If you see unusual MCMC fits, you can force the legacy (combination-based) inclusion-exclusion path: set env `MOIRE_USE_LEGACY_PAM=1` **before** starting R (e.g. `MOIRE_USE_LEGACY_PAM=1 Rscript -e "moire::run_mcmc(...)"` or set it in the shell before `R` / `devtools::load_all()`). Gray and combination implement the same formula; order of summation can cause small floating-point differences. The benchmark runs a Gray vs combination check (tol 1e-5) before timing.
- **MCMC C++ profiling (dev)**: To profile the MCMC hot path (e.g. inclusion-exclusion) on simulated data, build with `PKG_CXXFLAGS += -DMOIRE_ENABLE_PROFILER_REGISTRY` in `src/Makevars`, reinstall, then run `Rscript inst/scripts/profile_mcmc.R [minimal|small|vignette]`. Presets: `minimal` (tiny long-form data), `small` (20 samples × 10 loci), `vignette` (100 × 100). Env vars `PROFILE_BURNIN`, `PROFILE_SAMPLES`, `PROFILE_SEED` override defaults. Without the profiler flag, the script still reports wall-clock time.

#### MCMC profiling results (example: small preset)

Run: **small** (20 samples × 10 loci), 200 burnin + 200 samples, seed 42.

| Metric | Value |
|--------|--------|
| **Wall-clock** | ~19–34 s (elapsed; varies with build) |
| **Total C++ profiled** | ~104 s (sum of all scopes) |

**Top C++ hotspots (% of total profiled time):**

| Key | % | Calls | total_ms | avg_ms/call |
|-----|---|-------|----------|-------------|
| `Chain::calculate_transmission_likelihood` | 14–15 | 514k | ~15.6k–29.2k | 0.03–0.06 |
| `Chain::calc_transmission_process::all` | 13–14 | 514k | ~14.5k–27.3k | 0.03–0.05 |
| `Chain::calc_transmission_process::loop` | 12–12 | 514k | ~13.1k–25.1k | 0.025–0.049 |
| **`prob_any_missing::vectorized_cached::mobius_miss`** | **~11** | **~407k** | **~22.7k** | **~0.056** |
| `Chain::updates (burnin/sample)` | ~8–9 each | 200 each | ~9.5k–16.9k | 47–84 |
| `Chain::update_p` | ~7 | 400 | ~7.2k–13.7k | 18–34 |
| `Chain::update_samples` | ~3 | 400 | ~3.1k–5.4k | 7.7–13.4 |
| `Chain::update_m` / `update_r` / `update_eff_coi` | ~2–3 each | 400 | ~2.5k–5.1k | 6–13 |

**Per-algorithm breakdown (P(any missing) in transmission process):**

- **Möbius (cache miss)**: ~407k calls, **~22.7 s**, **~10.8%** of C++ time — dominant when `q` changes (e.g. new locus/sample).
- **Cache hit**: ~108k calls, ~54 ms total — negligible; cache is effective when same `q` is reused across samples.
- **Gray-code path**: not used in the run above (chain had been using cached Möbius).

**Recommendation:** Use **Gray code only** in the MCMC hot path. Benchmarks show Gray is ~8–10× faster per call than Möbius; in the profiled run ~79% of calls were cache misses (Möbius), so the cache did not pay off. The chain is configured to call `vectorized()` (Gray) only; the cache and Möbius path remain available in `prob_any_missing` for benchmarking or optional use.

#### Before/after: Gray-only vs Möbius+cache (small preset, 200 burnin + 200 samples)

| Metric | Möbius + cache (previous) | Gray only (current) | Speedup |
|--------|---------------------------|----------------------|---------|
| **Wall-clock (elapsed)** | ~19–33 s | **6.28 s** | **~3–5×** |
| **Chain::calc_transmission_process::loop** | 25,077 ms (0.049 ms/call) | **928 ms** (0.0018 ms/call) | **~27×** |
| **Chain::calculate_transmission_likelihood** | 29,174 ms | **2,622 ms** | **~11×** |
| **P(any missing) cost** | mobius_miss 22,740 ms + cache_hit 54 ms | included in loop above | — |

Gray-only removes the Möbius/cache overhead entirely; the transmission loop and likelihood are correspondingly much faster.

- **API compatibility**: The explicit `parallel_*` API (e.g. `parallel_reduce`, `parallel_sum`, `parallel_logsumexp`) is retained; no user-facing breaking changes in the current rollout.

#### Observation likelihood parallelization

Observation likelihood (`Chain::calculate_observation_likelihood` / `Chain::calc_observation_process`) is parallelized over (sample, locus) or locus where work is independent: in `initialize_likelihood` (parallel_for_2d over samples × loci), in `update_eps_pos`, `update_eps_neg`, `update_eff_coi`, `update_m`, `update_m_r`, and `update_samples` (parallel_for over loci for observation and for save/restore). Latent-genotype sampling remains sequential (RNG order). Parallelization uses the same thresholds as transmission (`MOIRE_PARALLEL_FOR_THRESHOLD`, `MOIRE_PARALLEL_FOR_2D_THRESHOLD`), so small datasets stay sequential and multi-chain runs respect `disable_nested_parallelism`. **Reproducibility:** parallel execution can change floating-point order; the same R seed may not give bit-identical MCMC output when observation runs in parallel.

## 🐛 **Debugging**

### **Debug Build**
```bash
CMAKE_BUILD_TYPE=Debug ./build_tests.sh
```

### **Debug Features**
- **Symbol Information**: Full debug symbols
- **Sanitizers**: Address sanitizer, undefined behavior sanitizer
- **Memory Leak Detection**: Automatic memory leak detection
- **Performance Profiling**: Detailed timing information

### **Coverage Analysis**
```bash
# Generate coverage report
TEST_TYPE=coverage ./run_tests.sh

# View coverage data
cd build
gcov *.gcno
lcov --capture --directory . --output-file coverage.info
genhtml coverage.info --output-directory coverage_html
```

## 🔄 **CI/CD Integration**

### **GitHub Actions**
```yaml
- name: Build C++ Tests
  run: |
    cd cpp_tests
    ./build_tests.sh
    
- name: Run C++ Tests
  run: |
    cd cpp_tests
    ./run_tests.sh
```

### **Local CI**
```bash
# Run full CI pipeline
just ci-full

# Run C++ specific CI
just cpp-build && just cpp-test
```

## 📝 **Adding New Tests**

### **1. Create Test File**
```cpp
#include "test_framework.hpp"
#include "your_header.h"

class YourTestSuite : public TestSuite {
public:
    YourTestSuite() : TestSuite("YourTestSuite") {
        add_test("Test Name", [this]() { test_function(); });
    }

private:
    void test_function() {
        // Your test implementation
    }
};
```

### **2. Add to CMakeLists.txt**
```cmake
add_executable(your_tests your_test.cpp)
target_link_libraries(your_tests 
    PRIVATE 
    Catch2::Catch2WithMain
    Threads::Threads
    TBB::tbb
    OpenMP::OpenMP_CXX
)
```

### **3. Update Build Scripts**
Add your test to the appropriate build and run scripts.

## 🎯 **Best Practices**

### **Test Design**
- **Single Responsibility**: Each test should test one specific behavior
- **Clear Naming**: Use descriptive test names
- **Independent Tests**: Tests should not depend on each other
- **Fast Execution**: Keep tests fast for quick feedback

### **Assertions**
- **Specific Assertions**: Use the most specific assertion macro
- **Clear Messages**: Provide clear failure messages
- **Edge Cases**: Test boundary conditions and error cases

### **Performance Testing**
- **Warmup Runs**: Include warmup runs to stabilize performance
- **Multiple Iterations**: Run multiple iterations for statistical significance
- **Realistic Data**: Use realistic test data sizes
- **Baseline Comparison**: Compare against known good baselines

### **Maintenance**
- **Regular Updates**: Keep tests updated with code changes
- **Documentation**: Document complex test logic
- **Review**: Regularly review and refactor tests
- **Coverage**: Maintain high test coverage

## 🚨 **Troubleshooting**

### **Common Issues**

1. **Build Failures**
   ```bash
   # Clean and rebuild
   rm -rf build
   ./build_tests.sh
   ```

2. **Test Failures**
   ```bash
   # Run with verbose output
   ./build/multivector_tests --verbose
   ```

3. **Performance Issues**
   ```bash
   # Check system resources
   htop
   # Run with fewer threads
   TBB_NUM_THREADS=2 ./run_tests.sh
   ```

4. **Memory Issues**
   ```bash
   # Run with memory debugging
   CMAKE_BUILD_TYPE=Debug ./build_tests.sh
   valgrind ./build/multivector_tests
   ```

### **Getting Help**
- Check build logs for compilation errors
- Use debug builds for detailed error information
- Run individual tests to isolate issues
- Check system resources and dependencies

This testing framework provides a robust foundation for C++ testing in the Moire package, ensuring code quality, performance, and reliability.

