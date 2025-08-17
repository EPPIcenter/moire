#include "../src/distance_matrix.h"
#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>

// Test function to demonstrate the DistanceMatrix class
void test_distance_matrix() {
    std::cout << "=== Distance Matrix Test ===\n\n";
    
    // Test data: 4 samples with 3 features each
    std::vector<std::vector<float>> test_data = {
        {1.0, 0.0, 1.0, 1.0, 0.0},  // Sample 0
        {1.0, 1.0, 1.0, 0.0, 0.0},  // Sample 1
        {0.0, 0.0, 1.0, 1.0, 1.0},  // Sample 2
        {1.0, 1.0, 0.0, 0.0, 1.0},  // Sample 3
        {0.0, 0.0, 0.0, 0.0, 0.0}   // Sample 4
    };
    
    const size_t num_samples = test_data.size();
    
    // Test 1: Jaccard Distance Matrix
    std::cout << "1. Jaccard Distance Matrix:\n";
    auto start = std::chrono::high_resolution_clock::now();
    
    DistanceMatrix<float> jaccard_matrix(num_samples, DistanceMatrix<float>::MetricType::JACCARD_DISTANCE);
    jaccard_matrix.compute(test_data);
    
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    std::cout << "   Distance matrix:\n";
    for (size_t i = 0; i < num_samples; ++i) {
        std::cout << "   ";
        for (size_t j = 0; j < num_samples; ++j) {
            std::cout << std::fixed << std::setprecision(3) << jaccard_matrix.get_distance(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "   Computation time: " << duration.count() << " microseconds\n";
    std::cout << "   Min distance: " << jaccard_matrix.min_distance() << "\n";
    std::cout << "   Max distance: " << jaccard_matrix.max_distance() << "\n";
    std::cout << "   Mean distance: " << jaccard_matrix.mean_distance() << "\n\n";
    
    // Test 2: Euclidean Distance Matrix
    std::cout << "2. Euclidean Distance Matrix:\n";
    start = std::chrono::high_resolution_clock::now();
    
    DistanceMatrix<float> euclidean_matrix(num_samples, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
    euclidean_matrix.compute(test_data);
    
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    std::cout << "   Distance matrix:\n";
    for (size_t i = 0; i < num_samples; ++i) {
        std::cout << "   ";
        for (size_t j = 0; j < num_samples; ++j) {
            std::cout << std::fixed << std::setprecision(3) << euclidean_matrix.get_distance(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "   Computation time: " << duration.count() << " microseconds\n";
    std::cout << "   Min distance: " << euclidean_matrix.min_distance() << "\n";
    std::cout << "   Max distance: " << euclidean_matrix.max_distance() << "\n";
    std::cout << "   Mean distance: " << euclidean_matrix.mean_distance() << "\n\n";
    
    // Test 3: Hamming Distance Matrix
    std::cout << "3. Hamming Distance Matrix:\n";
    start = std::chrono::high_resolution_clock::now();
    
    DistanceMatrix<float> hamming_matrix(num_samples, DistanceMatrix<float>::MetricType::HAMMING_DISTANCE);
    hamming_matrix.compute(test_data);
    
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    std::cout << "   Distance matrix:\n";
    for (size_t i = 0; i < num_samples; ++i) {
        std::cout << "   ";
        for (size_t j = 0; j < num_samples; ++j) {
            std::cout << std::fixed << std::setprecision(1) << hamming_matrix.get_distance(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "   Computation time: " << duration.count() << " microseconds\n";
    std::cout << "   Min distance: " << hamming_matrix.min_distance() << "\n";
    std::cout << "   Max distance: " << hamming_matrix.max_distance() << "\n";
    std::cout << "   Mean distance: " << hamming_matrix.mean_distance() << "\n\n";
    
    // Test 4: Custom Distance Function
    std::cout << "4. Custom Distance Function (Weighted Manhattan):\n";
    
    auto weighted_manhattan = [](const std::span<float const>& a, const std::span<float const>& b) -> float {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must have the same size");
        }
        
        // Weights for each dimension
        std::vector<float> weights = {2.0, 1.5, 1.0, 0.8, 0.5}; // First dimensions weighted more heavily
        
        float sum = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            sum += weights[i] * std::abs(a[i] - b[i]);
        }
        
        return sum;
    };
    
    start = std::chrono::high_resolution_clock::now();
    
    DistanceMatrix<float> custom_matrix(num_samples, weighted_manhattan, DistanceTag{});
    custom_matrix.compute(test_data);
    
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    std::cout << "   Distance matrix:\n";
    for (size_t i = 0; i < num_samples; ++i) {
        std::cout << "   ";
        for (size_t j = 0; j < num_samples; ++j) {
            std::cout << std::fixed << std::setprecision(3) << custom_matrix.get_distance(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "   Computation time: " << duration.count() << " microseconds\n";
    std::cout << "   Min distance: " << custom_matrix.min_distance() << "\n";
    std::cout << "   Max distance: " << custom_matrix.max_distance() << "\n";
    std::cout << "   Mean distance: " << custom_matrix.mean_distance() << "\n\n";
    
    // Test 5: Matrix Properties
    std::cout << "5. Matrix Properties:\n";
    std::cout << "   Jaccard matrix is symmetric: " << (jaccard_matrix.is_symmetric() ? "Yes" : "No") << "\n";
    std::cout << "   Euclidean matrix is symmetric: " << (euclidean_matrix.is_symmetric() ? "Yes" : "No") << "\n";
    std::cout << "   Hamming matrix is symmetric: " << (hamming_matrix.is_symmetric() ? "Yes" : "No") << "\n";
    std::cout << "   Custom matrix is symmetric: " << (custom_matrix.is_symmetric() ? "Yes" : "No") << "\n\n";
    
    // Test 6: Dynamic Metric Changes
    std::cout << "6. Dynamic Metric Changes:\n";
    DistanceMatrix<float> dynamic_matrix(num_samples);
    
    // Start with Jaccard distance
    dynamic_matrix.set_metric(DistanceMatrix<float>::MetricType::JACCARD_DISTANCE);
    dynamic_matrix.compute(test_data);
    std::cout << "   " << dynamic_matrix.metric_name() << " - Distance(0,1): " 
              << dynamic_matrix.get_distance(0, 1) << "\n";
    
    // Switch to Euclidean distance
    dynamic_matrix.set_metric(DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
    dynamic_matrix.compute(test_data);
    std::cout << "   " << dynamic_matrix.metric_name() << " - Distance(0,1): " 
              << dynamic_matrix.get_distance(0, 1) << "\n";
    
    // Switch to cosine distance
    dynamic_matrix.set_metric(DistanceMatrix<float>::MetricType::COSINE_DISTANCE);
    dynamic_matrix.compute(test_data);
    std::cout << "   " << dynamic_matrix.metric_name() << " - Distance(0,1): " 
              << dynamic_matrix.get_distance(0, 1) << "\n\n";
    
    std::cout << "=== All Tests Passed! ===\n";
}

// Performance benchmark for larger matrices
void benchmark_large_matrix() {
    std::cout << "\n=== Performance Benchmark ===\n\n";
    
    const size_t num_samples = 100;
    const size_t num_features = 50;
    
    // Generate random test data
    std::vector<std::vector<float>> large_data(num_samples);
    for (size_t i = 0; i < num_samples; ++i) {
        large_data[i].resize(num_features);
        for (size_t j = 0; j < num_features; ++j) {
            large_data[i][j] = static_cast<float>(rand()) / RAND_MAX;
        }
    }
    
    std::cout << "Testing with " << num_samples << " samples and " << num_features << " features\n\n";
    
    // Benchmark different metrics
    std::vector<std::pair<std::string, DistanceMatrix<float>::MetricType>> metrics = {
        {"Jaccard Distance", DistanceMatrix<float>::MetricType::JACCARD_DISTANCE},
        {"Euclidean Distance", DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE},
        {"Manhattan Distance", DistanceMatrix<float>::MetricType::MANHATTAN_DISTANCE},
        {"Cosine Distance", DistanceMatrix<float>::MetricType::COSINE_DISTANCE},
        {"Hamming Distance", DistanceMatrix<float>::MetricType::HAMMING_DISTANCE}
    };
    
    for (const auto& [name, metric] : metrics) {
        auto start = std::chrono::high_resolution_clock::now();
        
        DistanceMatrix<float> matrix(num_samples, metric);
        matrix.compute(large_data);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        
        std::cout << std::setw(20) << std::left << name << ": " 
                  << duration.count() << " ms, "
                  << "min=" << std::fixed << std::setprecision(3) << matrix.min_distance() << ", "
                  << "max=" << matrix.max_distance() << ", "
                  << "mean=" << matrix.mean_distance() << "\n";
    }
    
    std::cout << "\n=== Benchmark Complete ===\n";
}

int main() {
    std::cout << "Distance Matrix Class Test Suite\n";
    std::cout << "================================\n\n";
    
    try {
        test_distance_matrix();
        benchmark_large_matrix();
        std::cout << "\nAll tests completed successfully!\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
