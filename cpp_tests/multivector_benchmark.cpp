#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <numeric>
#include <execution>
#include <thread>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#include "multivector.h"


// Helper function to print benchmark results
void print_benchmark_result(const std::string& name, double time_ms) {
    std::cout << std::setw(30) << std::left << name << ": " 
              << std::fixed << std::setprecision(2) << time_ms << " ms" << std::endl;
}

// Helper function to run a benchmark
template<typename Func>
double run_benchmark(const std::string& name, Func func, int iterations = 5) {
    std::vector<double> times;
    
    // Multiple warm-up runs to stabilize the system
    for (int i = 0; i < 3; ++i) {
        func();
    }
    
    // Actual benchmark runs
    for (int i = 0; i < iterations; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
        times.push_back(duration);
    }
    
    // Calculate average time
    double avg_time = std::accumulate(times.begin(), times.end(), 0.0) / times.size();
    print_benchmark_result(name, avg_time);
    
    return avg_time;
}


// Test function to benchmark different reduction methods
void benchmark_reduce_functions() {
    std::cout << "\n=== MultiVector Reduce Benchmark ===\n" << std::endl;
    
    // Create a large MultiVector for testing
    const size_t dim1 = 100;
    const size_t dim2 = 100;
    const size_t dim3 = 100;
    
    std::array<size_t, 3> dimensions = {dim1, dim2, dim3};
    MultiVector<double, 3> mv(dimensions);
    
    // Fill with random data
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (size_t i = 0; i < dim1; ++i) {
        for (size_t j = 0; j < dim2; ++j) {
            for (size_t k = 0; k < dim3; ++k) {
                mv.at({i, j, k}) = dis(gen);
            }
        }
    }
    
    // Store results to prevent dead code elimination
    std::vector<double> results;
    results.reserve(10); // Reserve space for all results
    
    // Benchmark sequential reduction
    double seq_time = run_benchmark("Sequential reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Benchmark parallel reduction (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_time = run_benchmark("Parallel reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double par_unseq_time = run_benchmark("Parallel unsequenced reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par_unseq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double unseq_time = run_benchmark("Unsequenced reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::unseq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Benchmark built-in parallel implementation
    double builtin_par_time = run_benchmark("Built-in parallel reduce", [&]() {
        auto result = mv.parallel_reduce(std::plus<double>(), 0.0);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Verify results are consistent across execution policies
    std::cout << "\nResult verification:" << std::endl;
    bool results_consistent = true;
    for (size_t i = 1; i < results.size(); ++i) {
        if (std::abs(results[i] - results[0]) > 1e-10) {
            std::cout << "WARNING: Results are inconsistent! Sequential: " << results[0] 
                      << ", Other: " << results[i] << std::endl;
            results_consistent = false;
        }
    }
    if (results_consistent) {
        std::cout << "All results are consistent across execution policies." << std::endl;
    }
    
    // Calculate speedup
    double par_speedup = seq_time / par_time;
    double par_unseq_speedup = seq_time / par_unseq_time;
    double unseq_speedup = seq_time / unseq_time;
    double builtin_par_speedup = seq_time / builtin_par_time;
    
    std::cout << "\nSpeedup compared to sequential:" << std::endl;
    print_benchmark_result("Parallel speedup", par_speedup);
    print_benchmark_result("Parallel unsequenced speedup", par_unseq_speedup);
    print_benchmark_result("Unsequenced speedup", unseq_speedup);
    print_benchmark_result("Built-in parallel speedup", builtin_par_speedup);
#endif
    
    // Benchmark convenience methods
    double seq_sum_time = run_benchmark("Sequential sum", [&]() {
        auto result = mv.sum(std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_sum_time = run_benchmark("Parallel sum", [&]() {
        auto result = mv.sum(std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double builtin_par_sum_time = run_benchmark("Built-in parallel sum", [&]() {
        auto result = mv.parallel_sum();
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Verify sum results are consistent
    if (std::abs(results[results.size()-1] - results[results.size()-2]) > 1e-10) {
        std::cout << "WARNING: Sum results are inconsistent! Sequential: " << results[results.size()-2] 
                  << ", Parallel: " << results[results.size()-1] << std::endl;
    } else {
        std::cout << "Sum results are consistent across execution policies." << std::endl;
    }
    
    double sum_speedup = seq_sum_time / par_sum_time;
    double builtin_sum_speedup = seq_sum_time / builtin_par_sum_time;
    std::cout << "\nSum speedup:" << std::endl;
    print_benchmark_result("Parallel sum speedup", sum_speedup);
    print_benchmark_result("Built-in parallel sum speedup", builtin_sum_speedup);
#endif
}

// Test function to benchmark with different sizes
void benchmark_different_sizes() {
    std::cout << "\n=== MultiVector Size Benchmark ===\n" << std::endl;
    
    // Test with different sizes
    const std::vector<std::array<size_t, 3>> sizes = {
        {10, 10, 10},    // Small
        {50, 50, 50},    // Medium
        {100, 100, 100}, // Large
        {200, 200, 200}  // Very large
    };
    
    for (const auto& dims : sizes) {
        std::cout << "\nSize: " << dims[0] << "x" << dims[1] << "x" << dims[2] << std::endl;
        
        MultiVector<double, 3> mv(dims);
        
        // Fill with random data
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        for (size_t i = 0; i < dims[0]; ++i) {
            for (size_t j = 0; j < dims[1]; ++j) {
                for (size_t k = 0; k < dims[2]; ++k) {
                    mv.at({i, j, k}) = dis(gen);
                }
            }
        }
        
        // Store results to prevent dead code elimination
        std::vector<double> results;
        results.reserve(3); // Reserve space for results
        
        // Benchmark sequential reduction
        double seq_time = run_benchmark("Sequential reduce", [&]() {
            auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::seq);
            // Store a sample of the result to prevent optimization
            results.push_back(result.at({0, 0}));
        });
        
        // Benchmark parallel reduction (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
        double par_time = run_benchmark("Parallel reduce", [&]() {
            auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par);
            // Store a sample of the result to prevent optimization
            results.push_back(result.at({0, 0}));
        });
        
        // Benchmark built-in parallel implementation
        double builtin_par_time = run_benchmark("Built-in parallel reduce", [&]() {
            auto result = mv.parallel_reduce(std::plus<double>(), 0.0);
            // Store a sample of the result to prevent optimization
            results.push_back(result.at({0, 0}));
        });
        
        // Verify results are consistent
        bool results_consistent = true;
        for (size_t i = 1; i < results.size(); ++i) {
            if (std::abs(results[i] - results[0]) > 1e-10) {
                std::cout << "WARNING: Results are inconsistent! Sequential: " << results[0] 
                          << ", Other: " << results[i] << std::endl;
                results_consistent = false;
            }
        }
        if (results_consistent) {
            std::cout << "Results are consistent across execution policies." << std::endl;
        }
        
        // Calculate speedup
        double par_speedup = seq_time / par_time;
        double builtin_par_speedup = seq_time / builtin_par_time;
        std::cout << "Standard parallel speedup: " << std::fixed << std::setprecision(2) << par_speedup << "x" << std::endl;
        std::cout << "Built-in parallel speedup: " << std::fixed << std::setprecision(2) << builtin_par_speedup << "x" << std::endl;
#endif
    }
}

// Test function to benchmark different operations
void benchmark_different_operations() {
    std::cout << "\n=== MultiVector Operation Benchmark ===\n" << std::endl;
    
    // Create a MultiVector for testing
    const size_t dim1 = 100;
    const size_t dim2 = 100;
    const size_t dim3 = 100;
    
    std::array<size_t, 3> dimensions = {dim1, dim2, dim3};
    MultiVector<double, 3> mv(dimensions);
    
    // Fill with random data
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (size_t i = 0; i < dim1; ++i) {
        for (size_t j = 0; j < dim2; ++j) {
            for (size_t k = 0; k < dim3; ++k) {
                mv.at({i, j, k}) = dis(gen);
            }
        }
    }
    
    // Store results to prevent dead code elimination
    std::vector<double> results;
    results.reserve(12); // Reserve space for all results
    
    // Benchmark different operations with sequential execution
    double seq_sum_time = run_benchmark("Sequential sum", [&]() {
        auto result = mv.sum(std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double seq_product_time = run_benchmark("Sequential product", [&]() {
        auto result = mv.product(std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double seq_max_time = run_benchmark("Sequential max", [&]() {
        auto result = mv.max(std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double seq_min_time = run_benchmark("Sequential min", [&]() {
        auto result = mv.min(std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Benchmark different operations with parallel execution (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_sum_time = run_benchmark("Parallel sum", [&]() {
        auto result = mv.sum(std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double par_product_time = run_benchmark("Parallel product", [&]() {
        auto result = mv.product(std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double par_max_time = run_benchmark("Parallel max", [&]() {
        auto result = mv.max(std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double par_min_time = run_benchmark("Parallel min", [&]() {
        auto result = mv.min(std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Benchmark different operations with built-in parallel implementation
    double builtin_par_sum_time = run_benchmark("Built-in parallel sum", [&]() {
        auto result = mv.parallel_sum();
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double builtin_par_product_time = run_benchmark("Built-in parallel product", [&]() {
        auto result = mv.parallel_product();
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double builtin_par_max_time = run_benchmark("Built-in parallel max", [&]() {
        auto result = mv.parallel_max();
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double builtin_par_min_time = run_benchmark("Built-in parallel min", [&]() {
        auto result = mv.parallel_min();
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Verify results are consistent
    std::cout << "\nResult verification:" << std::endl;
    bool results_consistent = true;
    
    // Check sum results
    if (std::abs(results[4] - results[0]) > 1e-10) {
        std::cout << "WARNING: Sum results are inconsistent! Sequential: " << results[0] 
                  << ", Parallel: " << results[4] << std::endl;
        results_consistent = false;
    }
    
    // Check product results
    if (std::abs(results[5] - results[1]) > 1e-10) {
        std::cout << "WARNING: Product results are inconsistent! Sequential: " << results[1] 
                  << ", Parallel: " << results[5] << std::endl;
        results_consistent = false;
    }
    
    // Check max results
    if (std::abs(results[6] - results[2]) > 1e-10) {
        std::cout << "WARNING: Max results are inconsistent! Sequential: " << results[2] 
                  << ", Parallel: " << results[6] << std::endl;
        results_consistent = false;
    }
    
    // Check min results
    if (std::abs(results[7] - results[3]) > 1e-10) {
        std::cout << "WARNING: Min results are inconsistent! Sequential: " << results[3] 
                  << ", Parallel: " << results[7] << std::endl;
        results_consistent = false;
    }
    
    if (results_consistent) {
        std::cout << "All results are consistent across execution policies." << std::endl;
    }
    
    // Calculate speedups
    std::cout << "\nSpeedups compared to sequential:" << std::endl;
    print_benchmark_result("Sum speedup", seq_sum_time / par_sum_time);
    print_benchmark_result("Product speedup", seq_product_time / par_product_time);
    print_benchmark_result("Max speedup", seq_max_time / par_max_time);
    print_benchmark_result("Min speedup", seq_min_time / par_min_time);
    
    std::cout << "\nBuilt-in parallel speedups compared to sequential:" << std::endl;
    print_benchmark_result("Built-in sum speedup", seq_sum_time / builtin_par_sum_time);
    print_benchmark_result("Built-in product speedup", seq_product_time / builtin_par_product_time);
    print_benchmark_result("Built-in max speedup", seq_max_time / builtin_par_max_time);
    print_benchmark_result("Built-in min speedup", seq_min_time / builtin_par_min_time);
#endif
}

// Test function to benchmark 2D MultiVector
void benchmark_2d_multivector() {
    std::cout << "\n=== 2D MultiVector Benchmark ===\n" << std::endl;
    
    // Create a large 2D MultiVector for testing
    const size_t dim1 = 1000;
    const size_t dim2 = 1000;
    
    std::array<size_t, 2> dimensions = {dim1, dim2};
    MultiVector<double, 2> mv(dimensions);
    
    // Fill with random data
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (size_t i = 0; i < dim1; ++i) {
        for (size_t j = 0; j < dim2; ++j) {
            mv.at({i, j}) = dis(gen);
        }
    }
    
    // Store results to prevent dead code elimination
    std::vector<double> results;
    results.reserve(10); // Reserve space for all results
    
    // Benchmark sequential reduction
    double seq_time = run_benchmark("2D Sequential reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0}));
    });
    
    // Benchmark parallel reduction (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_time = run_benchmark("2D Parallel reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0}));
    });
    
    double par_unseq_time = run_benchmark("2D Parallel unsequenced reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par_unseq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0}));
    });
    
    double unseq_time = run_benchmark("2D Unsequenced reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::unseq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0}));
    });
    
    // Benchmark built-in parallel implementation
    double builtin_par_time = run_benchmark("2D Built-in parallel reduce", [&]() {
        auto result = mv.parallel_reduce(std::plus<double>(), 0.0);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0}));
    });
    
    // Verify results are consistent across execution policies
    std::cout << "\n2D Result verification:" << std::endl;
    bool results_consistent = true;
    for (size_t i = 1; i < results.size(); ++i) {
        // Use a more appropriate tolerance for floating-point comparison
        double relative_tolerance = 1e-10;
        double absolute_tolerance = 1e-10;
        double max_value = std::max(std::abs(results[0]), std::abs(results[i]));
        double tolerance = std::max(absolute_tolerance, relative_tolerance * max_value);
        
        if (std::abs(results[i] - results[0]) > tolerance) {
            std::cout << "WARNING: Results are inconsistent! Sequential: " << results[0] 
                      << ", Other: " << results[i] << " (diff: " << std::abs(results[i] - results[0]) << ")" << std::endl;
            results_consistent = false;
        }
    }
    if (results_consistent) {
        std::cout << "All 2D results are consistent across execution policies." << std::endl;
    }
    
    // Calculate speedup
    double par_speedup = seq_time / par_time;
    double par_unseq_speedup = seq_time / par_unseq_time;
    double unseq_speedup = seq_time / unseq_time;
    double builtin_par_speedup = seq_time / builtin_par_time;
    
    std::cout << "\n2D Speedup compared to sequential:" << std::endl;
    print_benchmark_result("2D Parallel speedup", par_speedup);
    print_benchmark_result("2D Parallel unsequenced speedup", par_unseq_speedup);
    print_benchmark_result("2D Unsequenced speedup", unseq_speedup);
    print_benchmark_result("2D Built-in parallel speedup", builtin_par_speedup);
#endif
    
    // Benchmark convenience methods
    double seq_sum_time = run_benchmark("2D Sequential sum", [&]() {
        auto result = mv.sum(std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0}));
    });
    
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_sum_time = run_benchmark("2D Parallel sum", [&]() {
        auto result = mv.sum(std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0}));
    });
    
    double builtin_par_sum_time = run_benchmark("2D Built-in parallel sum", [&]() {
        auto result = mv.parallel_sum();
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0}));
    });
    
    // Verify sum results are consistent
    double relative_tolerance = 1e-10;
    double absolute_tolerance = 1e-10;
    double max_value = std::max(std::abs(results[results.size()-2]), std::abs(results[results.size()-1]));
    double tolerance = std::max(absolute_tolerance, relative_tolerance * max_value);
    
    if (std::abs(results[results.size()-1] - results[results.size()-2]) > tolerance) {
        std::cout << "WARNING: 2D Sum results are inconsistent! Sequential: " << results[results.size()-2] 
                  << ", Parallel: " << results[results.size()-1] 
                  << " (diff: " << std::abs(results[results.size()-1] - results[results.size()-2]) << ")" << std::endl;
    } else {
        std::cout << "2D Sum results are consistent across execution policies." << std::endl;
    }
    
    double sum_speedup = seq_sum_time / par_sum_time;
    double builtin_sum_speedup = seq_sum_time / builtin_par_sum_time;
    std::cout << "\n2D Sum speedup:" << std::endl;
    print_benchmark_result("2D Parallel sum speedup", sum_speedup);
    print_benchmark_result("2D Built-in parallel sum speedup", builtin_sum_speedup);
#endif
}

// Test function to benchmark 4D MultiVector
void benchmark_4d_multivector() {
    std::cout << "\n=== 4D MultiVector Benchmark ===\n" << std::endl;
    
    // Create a 4D MultiVector for testing
    const size_t dim1 = 50;
    const size_t dim2 = 50;
    const size_t dim3 = 50;
    const size_t dim4 = 50;
    
    std::array<size_t, 4> dimensions = {dim1, dim2, dim3, dim4};
    MultiVector<double, 4> mv(dimensions);
    
    // Fill with random data
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (size_t i = 0; i < dim1; ++i) {
        for (size_t j = 0; j < dim2; ++j) {
            for (size_t k = 0; k < dim3; ++k) {
                for (size_t l = 0; l < dim4; ++l) {
                    mv.at({i, j, k, l}) = dis(gen);
                }
            }
        }
    }
    
    // Store results to prevent dead code elimination
    std::vector<double> results;
    results.reserve(10); // Reserve space for all results
    
    // Benchmark sequential reduction
    double seq_time = run_benchmark("4D Sequential reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Benchmark parallel reduction (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_time = run_benchmark("4D Parallel reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double par_unseq_time = run_benchmark("4D Parallel unsequenced reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par_unseq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double unseq_time = run_benchmark("4D Unsequenced reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::unseq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Benchmark built-in parallel implementation
    double builtin_par_time = run_benchmark("4D Built-in parallel reduce", [&]() {
        auto result = mv.parallel_reduce(std::plus<double>(), 0.0);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Verify results are consistent across execution policies
    std::cout << "\n4D Result verification:" << std::endl;
    bool results_consistent = true;
    for (size_t i = 1; i < results.size(); ++i) {
        // Use a more appropriate tolerance for floating-point comparison
        double relative_tolerance = 1e-10;
        double absolute_tolerance = 1e-10;
        double max_value = std::max(std::abs(results[0]), std::abs(results[i]));
        double tolerance = std::max(absolute_tolerance, relative_tolerance * max_value);
        
        if (std::abs(results[i] - results[0]) > tolerance) {
            std::cout << "WARNING: Results are inconsistent! Sequential: " << results[0] 
                      << ", Other: " << results[i] << " (diff: " << std::abs(results[i] - results[0]) << ")" << std::endl;
            results_consistent = false;
        }
    }
    if (results_consistent) {
        std::cout << "All 4D results are consistent across execution policies." << std::endl;
    }
    
    // Calculate speedup
    double par_speedup = seq_time / par_time;
    double par_unseq_speedup = seq_time / par_unseq_time;
    double unseq_speedup = seq_time / unseq_time;
    double builtin_par_speedup = seq_time / builtin_par_time;
    
    std::cout << "\n4D Speedup compared to sequential:" << std::endl;
    print_benchmark_result("4D Parallel speedup", par_speedup);
    print_benchmark_result("4D Parallel unsequenced speedup", par_unseq_speedup);
    print_benchmark_result("4D Unsequenced speedup", unseq_speedup);
    print_benchmark_result("4D Built-in parallel speedup", builtin_par_speedup);
#endif
    
    // Benchmark convenience methods
    double seq_sum_time = run_benchmark("4D Sequential sum", [&]() {
        auto result = mv.sum(std::execution::seq);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_sum_time = run_benchmark("4D Parallel sum", [&]() {
        auto result = mv.sum(std::execution::par);
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    double builtin_par_sum_time = run_benchmark("4D Built-in parallel sum", [&]() {
        auto result = mv.parallel_sum();
        // Store a sample of the result to prevent optimization
        results.push_back(result.at({0, 0}));
    });
    
    // Verify sum results are consistent
    double relative_tolerance = 1e-10;
    double absolute_tolerance = 1e-10;
    double max_value = std::max(std::abs(results[results.size()-2]), std::abs(results[results.size()-1]));
    double tolerance = std::max(absolute_tolerance, relative_tolerance * max_value);
    
    if (std::abs(results[results.size()-1] - results[results.size()-2]) > tolerance) {
        std::cout << "WARNING: 4D Sum results are inconsistent! Sequential: " << results[results.size()-2] 
                  << ", Parallel: " << results[results.size()-1] 
                  << " (diff: " << std::abs(results[results.size()-1] - results[results.size()-2]) << ")" << std::endl;
    } else {
        std::cout << "4D Sum results are consistent across execution policies." << std::endl;
    }
    
    double sum_speedup = seq_sum_time / par_sum_time;
    double builtin_sum_speedup = seq_sum_time / builtin_par_sum_time;
    std::cout << "\n4D Sum speedup:" << std::endl;
    print_benchmark_result("4D Parallel sum speedup", sum_speedup);
    print_benchmark_result("4D Built-in parallel sum speedup", builtin_sum_speedup);
#endif
}

// Test function to benchmark logsumexp operations
void benchmark_logsumexp() {
    std::cout << "\n=== MultiVector LogSumExp Benchmark ===\n" << std::endl;
    
    // Test with different sizes
    const std::vector<std::array<size_t, 3>> sizes = {
        {10, 10, 10},    // Small
        {50, 50, 50},    // Medium
        {100, 100, 100}, // Large
        {200, 200, 200}  // Very large
    };
    
    for (const auto& dims : sizes) {
        std::cout << "\nSize: " << dims[0] << "x" << dims[1] << "x" << dims[2] << std::endl;
        
        MultiVector<double, 3> mv(dims);
        
        // Fill with random data that won't cause overflow
        // Using values between -10 and 10 to avoid numerical issues
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(-10.0, 10.0);
        
        for (size_t i = 0; i < dims[0]; ++i) {
            for (size_t j = 0; j < dims[1]; ++j) {
                for (size_t k = 0; k < dims[2]; ++k) {
                    mv.at({i, j, k}) = dis(gen);
                }
            }
        }
        
        // Store results to prevent dead code elimination
        std::vector<double> results;
        results.reserve(3);
        
        // Benchmark sequential logsumexp
        double seq_time = run_benchmark("Sequential logsumexp", [&]() {
            auto result = mv.logsumexp(std::execution::seq);
            // Store a sample of the result to prevent optimization
            results.push_back(result.at({0, 0}));
        });
        
        // Benchmark parallel logsumexp
        double par_time = run_benchmark("Parallel logsumexp", [&]() {
            auto result = mv.parallel_logsumexp();
            // Store a sample of the result to prevent optimization
            results.push_back(result.at({0, 0}));
        });
        
        // Verify results are consistent
        bool results_consistent = true;
        for (size_t i = 1; i < results.size(); ++i) {
            if (std::abs(results[i] - results[0]) > 1e-10) {
                std::cout << "WARNING: Results are inconsistent! Sequential: " << results[0] 
                          << ", Parallel: " << results[i] << std::endl;
                results_consistent = false;
            }
        }
        if (results_consistent) {
            std::cout << "Results are consistent across execution policies." << std::endl;
        }
        
        // Calculate speedup
        double par_speedup = seq_time / par_time;
        std::cout << "Parallel speedup: " << std::fixed << std::setprecision(2) << par_speedup << "x" << std::endl;
        
        // Test numerical stability with challenging values
        MultiVector<double, 3> mv_challenge(std::array<size_t, 3>{2, 2, 4});
        
        // Fill with values that would cause overflow without the stability trick
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 2; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    mv_challenge.at({i, j, k}) = 100.0 + k;  // Values: 100, 101, 102, 103
                }
            }
        }
        
        std::cout << "\nTesting numerical stability..." << std::endl;
        
        // Benchmark stability test
        run_benchmark("Stability test (large values)", [&]() {
            auto result = mv_challenge.logsumexp();
            // The result should be approximately 103.44 for each slice
            // (as verified in the unit tests)
            for (size_t i = 0; i < 2; ++i) {
                for (size_t j = 0; j < 2; ++j) {
                    if (std::abs(result.at({i, j}) - 103.44) > 0.01) {
                        std::cout << "WARNING: Numerical stability issue detected!" << std::endl;
                    }
                }
            }
        });
    }
}

// Test function to benchmark 1D MultiVector
void benchmark_1d_multivector() {
    std::cout << "\n=== 1D MultiVector Benchmark ===\n" << std::endl;
    
    // Create a large 1D MultiVector for testing
    const size_t dim1 = 1000000; // Use a larger size for 1D to make the benchmark meaningful
    
    std::array<size_t, 1> dimensions = {dim1};
    MultiVector<double, 1> mv(dimensions);
    
    // Fill with random data
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (size_t i = 0; i < dim1; ++i) {
        mv.at({i}) = dis(gen);
    }
    
    // Store results to prevent dead code elimination
    std::vector<double> results;
    results.reserve(10); // Reserve space for all results
    
    // Benchmark sequential reduction
    double seq_time = run_benchmark("1D Sequential reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::seq);
        // Store the result to prevent optimization
        results.push_back(result);
    });
    
    // Benchmark parallel reduction (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_time = run_benchmark("1D Parallel reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par);
        // Store the result to prevent optimization
        results.push_back(result);
    });
    
    double par_unseq_time = run_benchmark("1D Parallel unsequenced reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::par_unseq);
        // Store the result to prevent optimization
        results.push_back(result);
    });
    
    double unseq_time = run_benchmark("1D Unsequenced reduce", [&]() {
        auto result = mv.reduce(std::plus<double>(), 0.0, std::execution::unseq);
        // Store the result to prevent optimization
        results.push_back(result);
    });
    
    // Benchmark built-in parallel implementation with TBB
    double builtin_par_time = run_benchmark("1D Built-in parallel reduce (TBB)", [&]() {
        auto result = mv.parallel_reduce(std::plus<double>(), 0.0);
        // Store the result to prevent optimization
        results.push_back(result);
    });
    
    // Verify results are consistent across execution policies
    std::cout << "\n1D Result verification:" << std::endl;
    bool results_consistent = true;
    for (size_t i = 1; i < results.size(); ++i) {
        // Use a more appropriate tolerance for floating-point comparison
        double relative_tolerance = 1e-10;
        double absolute_tolerance = 1e-10;
        double max_value = std::max(std::abs(results[0]), std::abs(results[i]));
        double tolerance = std::max(absolute_tolerance, relative_tolerance * max_value);
        
        if (std::abs(results[i] - results[0]) > tolerance) {
            std::cout << "WARNING: Results are inconsistent! Sequential: " << results[0] 
                      << ", Other: " << results[i] << " (diff: " << std::abs(results[i] - results[0]) << ")" << std::endl;
            results_consistent = false;
        }
    }
    if (results_consistent) {
        std::cout << "All 1D results are consistent across execution policies." << std::endl;
    }
    
    // Calculate speedup
    double par_speedup = seq_time / par_time;
    double par_unseq_speedup = seq_time / par_unseq_time;
    double unseq_speedup = seq_time / unseq_time;
    double builtin_par_speedup = seq_time / builtin_par_time;
    
    std::cout << "\n1D Speedup compared to sequential:" << std::endl;
    print_benchmark_result("1D Parallel speedup", par_speedup);
    print_benchmark_result("1D Parallel unsequenced speedup", par_unseq_speedup);
    print_benchmark_result("1D Unsequenced speedup", unseq_speedup);
    print_benchmark_result("1D Built-in parallel speedup (TBB)", builtin_par_speedup);
#endif
    
    // Benchmark convenience methods
    double seq_sum_time = run_benchmark("1D Sequential sum", [&]() {
        auto result = mv.sum(std::execution::seq);
        // Store the result to prevent optimization
        results.push_back(result);
    });
    
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_sum_time = run_benchmark("1D Parallel sum", [&]() {
        auto result = mv.sum(std::execution::par);
        // Store the result to prevent optimization
        results.push_back(result);
    });
    
    double builtin_par_sum_time = run_benchmark("1D Built-in parallel sum (TBB)", [&]() {
        auto result = mv.parallel_sum();
        // Store the result to prevent optimization
        results.push_back(result);
    });
    
    // Verify sum results are consistent
    double relative_tolerance = 1e-10;
    double absolute_tolerance = 1e-10;
    double max_value = std::max(std::abs(results[results.size()-2]), std::abs(results[results.size()-1]));
    double tolerance = std::max(absolute_tolerance, relative_tolerance * max_value);
    
    if (std::abs(results[results.size()-1] - results[results.size()-2]) > tolerance) {
        std::cout << "WARNING: 1D Sum results are inconsistent! Sequential: " << results[results.size()-2] 
                  << ", Parallel: " << results[results.size()-1] 
                  << " (diff: " << std::abs(results[results.size()-1] - results[results.size()-2]) << ")" << std::endl;
    } else {
        std::cout << "1D Sum results are consistent across execution policies." << std::endl;
    }
    
    double sum_speedup = seq_sum_time / par_sum_time;
    double builtin_sum_speedup = seq_sum_time / builtin_par_sum_time;
    std::cout << "\n1D Sum speedup:" << std::endl;
    print_benchmark_result("1D Parallel sum speedup", sum_speedup);
    print_benchmark_result("1D Built-in parallel sum speedup (TBB)", builtin_sum_speedup);
#endif
}

// Test function to benchmark element_transform operations
void benchmark_element_transform() {
    std::cout << "\n=== MultiVector Element Transform Benchmark ===\n" << std::endl;
    
    // Create a MultiVector for testing
    const size_t dim1 = 100;
    const size_t dim2 = 100;
    const size_t dim3 = 100;
    
    std::array<size_t, 3> dimensions = {dim1, dim2, dim3};
    MultiVector<double, 3> mv(dimensions);
    
    // Fill with random data
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (size_t i = 0; i < dim1; ++i) {
        for (size_t j = 0; j < dim2; ++j) {
            for (size_t k = 0; k < dim3; ++k) {
                mv.at({i, j, k}) = dis(gen);
            }
        }
    }
    
    // Create a copy for testing
    MultiVector<double, 3> mv_copy = mv;
    
    // Create values for each element position in the inner dimension
    std::vector<double> element_values(dim3);
    for (size_t i = 0; i < dim3; ++i) {
        element_values[i] = dis(gen);
    }
    
    // Store results to prevent dead code elimination
    std::vector<double> results;
    results.reserve(8); // Reserve space for all results
    
    // Benchmark sequential element_add
    double seq_element_add_time = run_benchmark("Sequential element_add", [&]() {
        mv.element_add(element_values);
        // Store a sample of the result to prevent optimization
        results.push_back(mv.at({0, 0, 0}));
        // Reset for next test
        mv = mv_copy;
    });
    
    // Benchmark sequential element_subtract
    double seq_element_subtract_time = run_benchmark("Sequential element_subtract", [&]() {
        mv.element_subtract(element_values);
        // Store a sample of the result to prevent optimization
        results.push_back(mv.at({0, 0, 0}));
        // Reset for next test
        mv = mv_copy;
    });
    
    // Benchmark sequential element_multiply
    double seq_element_multiply_time = run_benchmark("Sequential element_multiply", [&]() {
        mv.element_multiply(element_values);
        // Store a sample of the result to prevent optimization
        results.push_back(mv.at({0, 0, 0}));
        // Reset for next test
        mv = mv_copy;
    });
    
    // Benchmark sequential element_divide
    double seq_element_divide_time = run_benchmark("Sequential element_divide", [&]() {
        mv.element_divide(element_values);
        // Store a sample of the result to prevent optimization
        results.push_back(mv.at({0, 0, 0}));
        // Reset for next test
        mv = mv_copy;
    });
    
    // Benchmark parallel element operations (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    // Benchmark parallel_element_add
    double par_element_add_time = run_benchmark("Parallel element_add", [&]() {
        mv.parallel_element_add(element_values);
        // Store a sample of the result to prevent optimization
        results.push_back(mv.at({0, 0, 0}));
        // Reset for next test
        mv = mv_copy;
    });
    
    // Benchmark parallel_element_subtract
    double par_element_subtract_time = run_benchmark("Parallel element_subtract", [&]() {
        mv.parallel_element_subtract(element_values);
        // Store a sample of the result to prevent optimization
        results.push_back(mv.at({0, 0, 0}));
        // Reset for next test
        mv = mv_copy;
    });
    
    // Benchmark parallel_element_multiply
    double par_element_multiply_time = run_benchmark("Parallel element_multiply", [&]() {
        mv.parallel_element_multiply(element_values);
        // Store a sample of the result to prevent optimization
        results.push_back(mv.at({0, 0, 0}));
        // Reset for next test
        mv = mv_copy;
    });
    
    // Benchmark parallel_element_divide
    double par_element_divide_time = run_benchmark("Parallel element_divide", [&]() {
        mv.parallel_element_divide(element_values);
        // Store a sample of the result to prevent optimization
        results.push_back(mv.at({0, 0, 0}));
        // Reset for next test
        mv = mv_copy;
    });
    
    // Verify results are consistent
    std::cout << "\nResult verification:" << std::endl;
    bool results_consistent = true;
    
    // Check element_add results
    if (std::abs(results[4] - results[0]) > 1e-10) {
        std::cout << "WARNING: Element add results are inconsistent! Sequential: " << results[0] 
                  << ", Parallel: " << results[4] << std::endl;
        results_consistent = false;
    }
    
    // Check element_subtract results
    if (std::abs(results[5] - results[1]) > 1e-10) {
        std::cout << "WARNING: Element subtract results are inconsistent! Sequential: " << results[1] 
                  << ", Parallel: " << results[5] << std::endl;
        results_consistent = false;
    }
    
    // Check element_multiply results
    if (std::abs(results[6] - results[2]) > 1e-10) {
        std::cout << "WARNING: Element multiply results are inconsistent! Sequential: " << results[2] 
                  << ", Parallel: " << results[6] << std::endl;
        results_consistent = false;
    }
    
    // Check element_divide results
    if (std::abs(results[7] - results[3]) > 1e-10) {
        std::cout << "WARNING: Element divide results are inconsistent! Sequential: " << results[3] 
                  << ", Parallel: " << results[7] << std::endl;
        results_consistent = false;
    }
    
    if (results_consistent) {
        std::cout << "All element transform results are consistent across execution policies." << std::endl;
    }
    
    // Calculate speedups
    std::cout << "\nSpeedups compared to sequential:" << std::endl;
    print_benchmark_result("Element add speedup", seq_element_add_time / par_element_add_time);
    print_benchmark_result("Element subtract speedup", seq_element_subtract_time / par_element_subtract_time);
    print_benchmark_result("Element multiply speedup", seq_element_multiply_time / par_element_multiply_time);
    print_benchmark_result("Element divide speedup", seq_element_divide_time / par_element_divide_time);
#endif
    
    // Benchmark with different sizes
    std::cout << "\nBenchmarking element_transform with different sizes:" << std::endl;
    
    const std::vector<std::array<size_t, 3>> sizes = {
        {10, 10, 10},    // Small
        {50, 50, 50},    // Medium
        {100, 100, 100}, // Large
        {200, 200, 200}  // Very large
    };
    
    for (const auto& dims : sizes) {
        std::cout << "\nSize: " << dims[0] << "x" << dims[1] << "x" << dims[2] << std::endl;
        
        MultiVector<double, 3> mv_sized(dims);
        
        // Fill with random data
        for (size_t i = 0; i < dims[0]; ++i) {
            for (size_t j = 0; j < dims[1]; ++j) {
                for (size_t k = 0; k < dims[2]; ++k) {
                    mv_sized.at({i, j, k}) = dis(gen);
                }
            }
        }
        
        // Create a copy for testing
        MultiVector<double, 3> mv_sized_copy = mv_sized;
        
        // Create values for each element position in the inner dimension
        std::vector<double> sized_element_values(dims[2]);
        for (size_t i = 0; i < dims[2]; ++i) {
            sized_element_values[i] = dis(gen);
        }
        
        // Benchmark sequential element_add
        double sized_seq_element_add_time = run_benchmark("Sequential element_add", [&]() {
            mv_sized.element_add(sized_element_values);
            // Reset for next test
            mv_sized = mv_sized_copy;
        });
        
        // Benchmark parallel element_add (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
        double sized_par_element_add_time = run_benchmark("Parallel element_add", [&]() {
            mv_sized.parallel_element_add(sized_element_values);
            // Reset for next test
            mv_sized = mv_sized_copy;
        });
        
        // Calculate speedup
        double sized_speedup = sized_seq_element_add_time / sized_par_element_add_time;
        std::cout << "Element add speedup: " << std::fixed << std::setprecision(2) << sized_speedup << "x" << std::endl;
#endif
    }
    
    // Benchmark N=1 specialization
    std::cout << "\nBenchmarking N=1 specialization:" << std::endl;
    
    const size_t dim1d = 1000000; // Use a larger size for 1D to make the benchmark meaningful
    
    std::array<size_t, 1> dimensions1d = {dim1d};
    MultiVector<double, 1> mv1d(dimensions1d);
    
    // Fill with random data
    for (size_t i = 0; i < dim1d; ++i) {
        mv1d.at({i}) = dis(gen);
    }
    
    // Create a copy for testing
    MultiVector<double, 1> mv1d_copy = mv1d;
    
    // Create values for each element position
    std::vector<double> element_values1d(dim1d);
    for (size_t i = 0; i < dim1d; ++i) {
        element_values1d[i] = dis(gen);
    }
    
    // Benchmark sequential element_add
    double seq_element_add_time_1d = run_benchmark("1D Sequential element_add", [&]() {
        mv1d.element_add(element_values1d);
        // Reset for next test
        mv1d = mv1d_copy;
    });
    
    // Benchmark parallel element_add (if available)
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    double par_element_add_time_1d = run_benchmark("1D Parallel element_add", [&]() {
        mv1d.parallel_element_add(element_values1d);
        // Reset for next test
        mv1d = mv1d_copy;
    });
    
    // Calculate speedup
    double speedup_1d = seq_element_add_time_1d / par_element_add_time_1d;
    std::cout << "1D Element add speedup: " << std::fixed << std::setprecision(2) << speedup_1d << "x" << std::endl;
#endif
}

// Test function to benchmark elementwise operations between multivectors
void benchmark_elementwise_operations() {
    std::cout << "\n=== MultiVector Elementwise Operations Benchmark ===\n" << std::endl;
    
    // Create two MultiVectors for testing
    const size_t dim1 = 100;
    const size_t dim2 = 100;
    const size_t dim3 = 100;
    
    std::array<size_t, 3> dimensions = {dim1, dim2, dim3};
    MultiVector<double, 3> mv1(dimensions);
    MultiVector<double, 3> mv2(dimensions);
    
    // Fill with random data
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    
    for (size_t i = 0; i < dim1; ++i) {
        for (size_t j = 0; j < dim2; ++j) {
            for (size_t k = 0; k < dim3; ++k) {
                mv1.at({i, j, k}) = dis(gen);
                mv2.at({i, j, k}) = dis(gen);
            }
        }
    }
    
    // Store results to prevent dead code elimination
    std::vector<double> results;
    results.reserve(8); // Reserve space for all results
    
    // Benchmark addition
    double add_time = run_benchmark("Elementwise addition", [&]() {
        MultiVector<double, 3> result = mv1 + mv2;
        results.push_back(result.at({0, 0, 0})); // Prevent dead code elimination
    });
    
    // Benchmark subtraction
    double subtract_time = run_benchmark("Elementwise subtraction", [&]() {
        MultiVector<double, 3> result = mv2 - mv1;
        results.push_back(result.at({0, 0, 0})); // Prevent dead code elimination
    });
    
    // Benchmark multiplication
    double multiply_time = run_benchmark("Elementwise multiplication", [&]() {
        MultiVector<double, 3> result = mv1 * mv2;
        results.push_back(result.at({0, 0, 0})); // Prevent dead code elimination
    });
    
    // Benchmark division
    double divide_time = run_benchmark("Elementwise division", [&]() {
        MultiVector<double, 3> result = mv2 / mv1;
        results.push_back(result.at({0, 0, 0})); // Prevent dead code elimination
    });
    
    // Benchmark in-place addition
    double inplace_add_time = run_benchmark("In-place addition", [&]() {
        MultiVector<double, 3> mv_copy = mv1;
        mv_copy += mv2;
        results.push_back(mv_copy.at({0, 0, 0})); // Prevent dead code elimination
    });
    
    // Benchmark in-place subtraction
    double inplace_subtract_time = run_benchmark("In-place subtraction", [&]() {
        MultiVector<double, 3> mv_copy = mv1;
        mv_copy -= mv2;
        results.push_back(mv_copy.at({0, 0, 0})); // Prevent dead code elimination
    });
    
    // Benchmark in-place multiplication
    double inplace_multiply_time = run_benchmark("In-place multiplication", [&]() {
        MultiVector<double, 3> mv_copy = mv1;
        mv_copy *= mv2;
        results.push_back(mv_copy.at({0, 0, 0})); // Prevent dead code elimination
    });
    
    // Benchmark in-place division
    double inplace_divide_time = run_benchmark("In-place division", [&]() {
        MultiVector<double, 3> mv_copy = mv1;
        mv_copy /= mv2;
        results.push_back(mv_copy.at({0, 0, 0})); // Prevent dead code elimination
    });
    
    // Print summary
    std::cout << "\nElementwise Operations Summary:" << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << "Addition: " << add_time << " ms" << std::endl;
    std::cout << "Subtraction: " << subtract_time << " ms" << std::endl;
    std::cout << "Multiplication: " << multiply_time << " ms" << std::endl;
    std::cout << "Division: " << divide_time << " ms" << std::endl;
    std::cout << "In-place Addition: " << inplace_add_time << " ms" << std::endl;
    std::cout << "In-place Subtraction: " << inplace_subtract_time << " ms" << std::endl;
    std::cout << "In-place Multiplication: " << inplace_multiply_time << " ms" << std::endl;
    std::cout << "In-place Division: " << inplace_divide_time << " ms" << std::endl;
    
    // Compare in-place vs non-in-place operations
    std::cout << "\nIn-place vs Non-in-place Comparison:" << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Addition: " << (add_time / inplace_add_time) << "x faster in-place" << std::endl;
    std::cout << "Subtraction: " << (subtract_time / inplace_subtract_time) << "x faster in-place" << std::endl;
    std::cout << "Multiplication: " << (multiply_time / inplace_multiply_time) << "x faster in-place" << std::endl;
    std::cout << "Division: " << (divide_time / inplace_divide_time) << "x faster in-place" << std::endl;
    
    // Benchmark with different sizes
    std::cout << "\nBenchmarking with different sizes:" << std::endl;
    std::cout << "----------------------------------" << std::endl;
    
    // Small size
    std::array<size_t, 3> small_dims = {10, 10, 10};
    MultiVector<double, 3> small_mv1(small_dims);
    MultiVector<double, 3> small_mv2(small_dims);
    
    for (size_t i = 0; i < small_dims[0]; ++i) {
        for (size_t j = 0; j < small_dims[1]; ++j) {
            for (size_t k = 0; k < small_dims[2]; ++k) {
                small_mv1.at({i, j, k}) = dis(gen);
                small_mv2.at({i, j, k}) = dis(gen);
            }
        }
    }
    
    double small_add_time = run_benchmark("Small size addition (10x10x10)", [&]() {
        MultiVector<double, 3> result = small_mv1 + small_mv2;
        results.push_back(result.at({0, 0, 0}));
    });
    
    // Medium size
    std::array<size_t, 3> medium_dims = {50, 50, 50};
    MultiVector<double, 3> medium_mv1(medium_dims);
    MultiVector<double, 3> medium_mv2(medium_dims);
    
    for (size_t i = 0; i < medium_dims[0]; ++i) {
        for (size_t j = 0; j < medium_dims[1]; ++j) {
            for (size_t k = 0; k < medium_dims[2]; ++k) {
                medium_mv1.at({i, j, k}) = dis(gen);
                medium_mv2.at({i, j, k}) = dis(gen);
            }
        }
    }
    
    double medium_add_time = run_benchmark("Medium size addition (50x50x50)", [&]() {
        MultiVector<double, 3> result = medium_mv1 + medium_mv2;
        results.push_back(result.at({0, 0, 0}));
    });
    
    // Large size
    std::array<size_t, 3> large_dims = {200, 200, 200};
    MultiVector<double, 3> large_mv1(large_dims);
    MultiVector<double, 3> large_mv2(large_dims);
    
    for (size_t i = 0; i < large_dims[0]; ++i) {
        for (size_t j = 0; j < large_dims[1]; ++j) {
            for (size_t k = 0; k < large_dims[2]; ++k) {
                large_mv1.at({i, j, k}) = dis(gen);
                large_mv2.at({i, j, k}) = dis(gen);
            }
        }
    }
    
    double large_add_time = run_benchmark("Large size addition (200x200x200)", [&]() {
        MultiVector<double, 3> result = large_mv1 + large_mv2;
        results.push_back(result.at({0, 0, 0}));
    });
    
    // Print size comparison
    std::cout << "\nSize Comparison:" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout << "Small (10x10x10): " << small_add_time << " ms" << std::endl;
    std::cout << "Medium (50x50x50): " << medium_add_time << " ms" << std::endl;
    std::cout << "Large (200x200x200): " << large_add_time << " ms" << std::endl;
    
    // Calculate scaling factors
    double small_to_medium = medium_add_time / small_add_time;
    double medium_to_large = large_add_time / medium_add_time;
    double small_to_large = large_add_time / small_add_time;
    
    std::cout << "\nScaling Factors:" << std::endl;
    std::cout << "----------------" << std::endl;
    std::cout << "Small to Medium: " << small_to_medium << "x (expected ~125x)" << std::endl;
    std::cout << "Medium to Large: " << medium_to_large << "x (expected ~64x)" << std::endl;
    std::cout << "Small to Large: " << small_to_large << "x (expected ~8000x)" << std::endl;
}

void benchmark_chained_reductions() {
    std::cout << "\n=== MultiVector Chained Reductions Benchmark ===\n" << std::endl;
    
    // Create a large 3D MultiVector
    std::array<size_t, 3> dimensions = {100, 100, 100};
    MultiVector<double, 3> mv(dimensions);
    
    // Fill with random data
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-10.0, 10.0);
    
    for (size_t i = 0; i < dimensions[0]; ++i) {
        for (size_t j = 0; j < dimensions[1]; ++j) {
            for (size_t k = 0; k < dimensions[2]; ++k) {
                mv.at({i, j, k}) = dis(gen);
            }
        }
    }
    
    // Benchmark chained operations
    run_benchmark("Chained sum->logsumexp->sum", [&]() {
        return mv.parallel_sum().parallel_logsumexp().parallel_sum();
    });
    
    // Benchmark combined approach (single pass)
    run_benchmark("Combined single-pass reduction", [&]() {
        // Create a 2D MultiVector to store intermediate results
        std::array<size_t, 2> dim2d = {dimensions[0], dimensions[1]};
        MultiVector<double, 2> sum1(dim2d);
        
        // First pass: compute sum along the innermost dimension
        for (size_t i = 0; i < dimensions[0]; ++i) {
            for (size_t j = 0; j < dimensions[1]; ++j) {
                double sum = 0.0;
                for (size_t k = 0; k < dimensions[2]; ++k) {
                    sum += mv.at({i, j, k});
                }
                sum1.at({i, j}) = sum;
            }
        }
        
        // Second pass: compute logsumexp along the innermost dimension
        std::array<size_t, 1> dim1d = {dimensions[0]};
        MultiVector<double, 1> logsumexp(dim1d);
        
        for (size_t i = 0; i < dimensions[0]; ++i) {
            double max_val = std::numeric_limits<double>::lowest();
            double sum = 0.0;
            
            // Find max for numerical stability
            for (size_t j = 0; j < dimensions[1]; ++j) {
                max_val = std::max(max_val, sum1.at({i, j}));
            }
            
            // Compute exp(x - max) and sum
            for (size_t j = 0; j < dimensions[1]; ++j) {
                sum += std::exp(sum1.at({i, j}) - max_val);
            }
            
            // Compute log(sum) + max
            logsumexp.at({i}) = std::log(sum) + max_val;
        }
        
        // Final pass: compute sum along the innermost dimension
        double sum2 = 0.0;
        for (size_t i = 0; i < dimensions[0]; ++i) {
            sum2 += logsumexp.at({i});
        }
        
        return sum2;
    });
}

int main() {
    std::cout << "MultiVector Benchmark Test" << std::endl;
    std::cout << "=========================" << std::endl;
    
    // Check if parallel execution is available
#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
    std::cout << "Parallel execution is available." << std::endl;
#else
    std::cout << "Parallel execution is not available." << std::endl;
#endif
    std::cout << std::endl;
    
    benchmark_1d_multivector();
    benchmark_2d_multivector();
    benchmark_reduce_functions();
    benchmark_4d_multivector();
    benchmark_different_sizes();
    benchmark_different_operations();
    benchmark_logsumexp();
    benchmark_element_transform();
    benchmark_elementwise_operations();
    benchmark_chained_reductions();
    
    std::cout << "\nBenchmark completed." << std::endl;
    return 0;
} 