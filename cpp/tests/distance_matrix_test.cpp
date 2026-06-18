#include "test_framework.hpp"
#include "distance_matrix.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <span>

namespace moire_test {

// Test suite for DistanceMatrix
class DistanceMatrixTestSuite : public TestSuite {
public:
    DistanceMatrixTestSuite() : TestSuite("DistanceMatrix") {
        add_test("Construction", [this]() { test_construction(); });
        add_test("Jaccard Distance", [this]() { test_jaccard_distance(); });
        add_test("Euclidean Distance", [this]() { test_euclidean_distance(); });
        add_test("Manhattan Distance", [this]() { test_manhattan_distance(); });
        add_test("Cosine Distance", [this]() { test_cosine_distance(); });
        add_test("Distance Properties", [this]() { test_distance_properties(); });
        add_test("Statistics", [this]() { test_statistics(); });
        add_test("Edge Cases", [this]() { test_edge_cases(); });
        add_test("Performance", [this]() { test_performance(); });
        
        // Additional tests from old framework
        add_test("Hamming Distance", [this]() { test_hamming_distance(); });
        add_test("Custom Distance Function", [this]() { test_custom_distance_function(); });
        add_test("Matrix Properties", [this]() { test_matrix_properties(); });
        add_test("Dynamic Metric Changes", [this]() { test_dynamic_metric_changes(); });
        add_test("Large Matrix Benchmark", [this]() { test_large_matrix_benchmark(); });
    }

private:
    void test_construction() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 0.0, 1.0},
            {0.0, 1.0, 0.0},
            {1.0, 1.0, 1.0}
        };
        
        DistanceMatrix<float> jaccard_matrix(3, DistanceMatrix<float>::MetricType::JACCARD_DISTANCE);
        DistanceMatrix<float> euclidean_matrix(3, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        
        ASSERT_EQ(jaccard_matrix.num_samples(), 3);
        ASSERT_EQ(euclidean_matrix.num_samples(), 3);
    }

    void test_jaccard_distance() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 0.0, 1.0, 1.0, 0.0},  // Sample 0
            {1.0, 1.0, 1.0, 0.0, 0.0},  // Sample 1
            {0.0, 0.0, 1.0, 1.0, 1.0},  // Sample 2
            {1.0, 1.0, 0.0, 0.0, 1.0},  // Sample 3
            {0.0, 0.0, 0.0, 0.0, 0.0}   // Sample 4
        };
        
        DistanceMatrix<float> matrix(5, DistanceMatrix<float>::MetricType::JACCARD_DISTANCE);
        matrix.compute(test_data);
        
        // Test diagonal (should be 0)
        for (size_t i = 0; i < 5; ++i) {
            ASSERT_EQ(matrix.get_distance(i, i), 0.0f);
        }
        
        // Test symmetry
        for (size_t i = 0; i < 5; ++i) {
            for (size_t j = 0; j < 5; ++j) {
                ASSERT_EQ(matrix.get_distance(i, j), matrix.get_distance(j, i));
            }
        }
        
        // Test specific known distances
        // Jaccard distance between identical vectors should be 0
        ASSERT_EQ(matrix.get_distance(0, 0), 0.0f);
        
        // Test non-negative distances
        for (size_t i = 0; i < 5; ++i) {
            for (size_t j = 0; j < 5; ++j) {
                ASSERT_GE(matrix.get_distance(i, j), 0.0f);
            }
        }
    }

    void test_euclidean_distance() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
        };
        
        DistanceMatrix<float> matrix(3, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        matrix.compute(test_data);
        
        // Test diagonal (should be 0)
        for (size_t i = 0; i < 3; ++i) {
            ASSERT_EQ(matrix.get_distance(i, i), 0.0f);
        }
        
        // Test symmetry
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                ASSERT_EQ(matrix.get_distance(i, j), matrix.get_distance(j, i));
            }
        }
        
        // Test specific known distances
        // Distance between (1,2,3) and (4,5,6) should be sqrt(27) ≈ 5.196
        float expected_distance = std::sqrt(27.0f);
        float actual_distance = matrix.get_distance(0, 1);
        ASSERT_LT(std::abs(actual_distance - expected_distance), 0.001f);
    }

    void test_manhattan_distance() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
        };
        
        DistanceMatrix<float> matrix(3, DistanceMatrix<float>::MetricType::MANHATTAN_DISTANCE);
        matrix.compute(test_data);
        
        // Test diagonal (should be 0)
        for (size_t i = 0; i < 3; ++i) {
            ASSERT_EQ(matrix.get_distance(i, i), 0.0f);
        }
        
        // Test symmetry
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                ASSERT_EQ(matrix.get_distance(i, j), matrix.get_distance(j, i));
            }
        }
        
        // Test specific known distances
        // Manhattan distance between (1,2,3) and (4,5,6) should be 9
        float expected_distance = 9.0f;
        float actual_distance = matrix.get_distance(0, 1);
        ASSERT_EQ(actual_distance, expected_distance);
    }

    void test_cosine_distance() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0},
            {1.0, 1.0, 1.0}
        };
        
        DistanceMatrix<float> matrix(4, DistanceMatrix<float>::MetricType::COSINE_DISTANCE);
        matrix.compute(test_data);
        
        // Test diagonal (should be 0)
        for (size_t i = 0; i < 4; ++i) {
            ASSERT_EQ(matrix.get_distance(i, i), 0.0f);
        }
        
        // Test symmetry
        for (size_t i = 0; i < 4; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                ASSERT_EQ(matrix.get_distance(i, j), matrix.get_distance(j, i));
            }
        }
        
        // Test orthogonality (cosine distance between orthogonal vectors should be 1)
        float orthogonal_distance = matrix.get_distance(0, 1);
        ASSERT_LT(std::abs(orthogonal_distance - 1.0f), 0.001f);
    }

    void test_distance_properties() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
        };
        
        DistanceMatrix<float> matrix(3, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        matrix.compute(test_data);
        
        // Test triangle inequality
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 3; ++k) {
                    float d_ij = matrix.get_distance(i, j);
                    float d_ik = matrix.get_distance(i, k);
                    float d_kj = matrix.get_distance(k, j);
                    ASSERT_LE(d_ij, d_ik + d_kj + 0.001f); // Small tolerance for floating point
                }
            }
        }
    }

    void test_statistics() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
        };
        
        DistanceMatrix<float> matrix(3, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        matrix.compute(test_data);
        
        // Test statistics
        float min_dist = matrix.min_distance();
        float max_dist = matrix.max_distance();
        float mean_dist = matrix.mean_distance();
        
        ASSERT_GE(min_dist, 0.0f);
        ASSERT_GE(max_dist, min_dist);
        ASSERT_GE(mean_dist, min_dist);
        ASSERT_LE(mean_dist, max_dist);
        
        // Note: This implementation computes actual minimum distance across all pairs
        // rather than setting diagonal elements to 0. The minimum distance is the
        // distance between the closest pair of samples, not necessarily 0.
    }

    void test_edge_cases() {
        // Test with single sample
        std::vector<std::vector<float>> single_sample = {{1.0, 2.0, 3.0}};
        DistanceMatrix<float> single_matrix(1, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        single_matrix.compute(single_sample);
        ASSERT_EQ(single_matrix.get_distance(0, 0), 0.0f);
        
        // Test with identical samples
        std::vector<std::vector<float>> identical_samples = {
            {1.0, 2.0, 3.0},
            {1.0, 2.0, 3.0},
            {1.0, 2.0, 3.0}
        };
        DistanceMatrix<float> identical_matrix(3, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        identical_matrix.compute(identical_samples);
        
        // All distances should be 0
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                ASSERT_EQ(identical_matrix.get_distance(i, j), 0.0f);
            }
        }
        
        // Test with zero vectors
        std::vector<std::vector<float>> zero_vectors = {
            {0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0}
        };
        DistanceMatrix<float> zero_matrix(2, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        zero_matrix.compute(zero_vectors);
        ASSERT_EQ(zero_matrix.get_distance(0, 1), 0.0f);
    }

    void test_performance() {
        // Create larger test data
        std::vector<std::vector<float>> large_data;
        for (size_t i = 0; i < 100; ++i) {
            std::vector<float> sample;
            for (size_t j = 0; j < 50; ++j) {
                sample.push_back(static_cast<float>(i + j));
            }
            large_data.push_back(sample);
        }
        
        DistanceMatrix<float> matrix(100, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        
        auto start = std::chrono::high_resolution_clock::now();
        matrix.compute(large_data);
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        // Performance should be reasonable (less than 5 seconds for this test)
        ASSERT_LT(duration, 5000);
        
        // Verify the computation worked
        ASSERT_EQ(matrix.num_samples(), 100);
        for (size_t i = 0; i < 100; ++i) {
            ASSERT_EQ(matrix.get_distance(i, i), 0.0f);
        }
    }

    void test_hamming_distance() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 0.0, 1.0, 1.0, 0.0},  // Sample 0
            {1.0, 1.0, 1.0, 0.0, 0.0},  // Sample 1
            {0.0, 0.0, 1.0, 1.0, 1.0},  // Sample 2
            {1.0, 1.0, 0.0, 0.0, 1.0},  // Sample 3
            {0.0, 0.0, 0.0, 0.0, 0.0}   // Sample 4
        };
        
        DistanceMatrix<float> matrix(5, DistanceMatrix<float>::MetricType::HAMMING_DISTANCE);
        matrix.compute(test_data);
        
        // Test diagonal (should be 0)
        for (size_t i = 0; i < 5; ++i) {
            ASSERT_EQ(matrix.get_distance(i, i), 0.0f);
        }
        
        // Test symmetry
        for (size_t i = 0; i < 5; ++i) {
            for (size_t j = 0; j < 5; ++j) {
                ASSERT_EQ(matrix.get_distance(i, j), matrix.get_distance(j, i));
            }
        }
        
        // Test specific known distances
        // Hamming distance between identical vectors should be 0
        ASSERT_EQ(matrix.get_distance(0, 0), 0.0f);
        
        // Test non-negative distances
        for (size_t i = 0; i < 5; ++i) {
            for (size_t j = 0; j < 5; ++j) {
                ASSERT_GE(matrix.get_distance(i, j), 0.0f);
            }
        }
    }

    void test_custom_distance_function() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 2.0, 3.0, 4.0, 5.0},
            {2.0, 3.0, 4.0, 5.0, 6.0},
            {3.0, 4.0, 5.0, 6.0, 7.0}
        };
        
        // Custom weighted Manhattan distance function
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
        
        DistanceMatrix<float> matrix(3, weighted_manhattan, DistanceTag{});
        matrix.compute(test_data);
        
        // Test diagonal (should be 0)
        for (size_t i = 0; i < 3; ++i) {
            ASSERT_EQ(matrix.get_distance(i, i), 0.0f);
        }
        
        // Test symmetry
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                ASSERT_EQ(matrix.get_distance(i, j), matrix.get_distance(j, i));
            }
        }
        
        // Test non-negative distances
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                ASSERT_GE(matrix.get_distance(i, j), 0.0f);
            }
        }
    }

    void test_matrix_properties() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 2.0, 3.0},
            {4.0, 5.0, 6.0},
            {7.0, 8.0, 9.0}
        };
        
        DistanceMatrix<float> jaccard_matrix(3, DistanceMatrix<float>::MetricType::JACCARD_DISTANCE);
        DistanceMatrix<float> euclidean_matrix(3, DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        DistanceMatrix<float> hamming_matrix(3, DistanceMatrix<float>::MetricType::HAMMING_DISTANCE);
        
        jaccard_matrix.compute(test_data);
        euclidean_matrix.compute(test_data);
        hamming_matrix.compute(test_data);
        
        // Test matrix properties
        ASSERT_TRUE(jaccard_matrix.is_symmetric());
        ASSERT_TRUE(euclidean_matrix.is_symmetric());
        ASSERT_TRUE(hamming_matrix.is_symmetric());
        
        // Test that diagonal elements are 0
        for (size_t i = 0; i < 3; ++i) {
            ASSERT_EQ(jaccard_matrix.get_distance(i, i), 0.0f);
            ASSERT_EQ(euclidean_matrix.get_distance(i, i), 0.0f);
            ASSERT_EQ(hamming_matrix.get_distance(i, i), 0.0f);
        }
    }

    void test_dynamic_metric_changes() {
        std::vector<std::vector<float>> test_data = {
            {1.0, 0.0, 1.0},
            {0.0, 1.0, 0.0},
            {1.0, 1.0, 1.0}
        };
        
        DistanceMatrix<float> dynamic_matrix(3);
        
        // Start with Jaccard distance
        dynamic_matrix.set_metric(DistanceMatrix<float>::MetricType::JACCARD_DISTANCE);
        dynamic_matrix.compute(test_data);
        float jaccard_distance = dynamic_matrix.get_distance(0, 1);
        ASSERT_GE(jaccard_distance, 0.0f);
        
        // Switch to Euclidean distance
        dynamic_matrix.set_metric(DistanceMatrix<float>::MetricType::EUCLIDEAN_DISTANCE);
        dynamic_matrix.compute(test_data);
        float euclidean_distance = dynamic_matrix.get_distance(0, 1);
        ASSERT_GE(euclidean_distance, 0.0f);
        
        // Switch to cosine distance
        dynamic_matrix.set_metric(DistanceMatrix<float>::MetricType::COSINE_DISTANCE);
        dynamic_matrix.compute(test_data);
        float cosine_distance = dynamic_matrix.get_distance(0, 1);
        ASSERT_GE(cosine_distance, 0.0f);
        
        // Verify that different metrics give different results
        ASSERT_NE(jaccard_distance, euclidean_distance);
        ASSERT_NE(euclidean_distance, cosine_distance);
    }

    void test_large_matrix_benchmark() {
        const size_t num_samples = 50;  // Reduced for test performance
        const size_t num_features = 20;
        
        // Generate test data
        std::vector<std::vector<float>> large_data(num_samples);
        for (size_t i = 0; i < num_samples; ++i) {
            large_data[i].resize(num_features);
            for (size_t j = 0; j < num_features; ++j) {
                large_data[i][j] = static_cast<float>(i + j) / 10.0f;
            }
        }
        
        // Test different metrics
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
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            
            // Performance should be reasonable (less than 2 seconds for this test)
            ASSERT_LT(duration, 2000);
            
            // Verify the computation worked
            ASSERT_EQ(matrix.num_samples(), num_samples);
            for (size_t i = 0; i < num_samples; ++i) {
                ASSERT_EQ(matrix.get_distance(i, i), 0.0f);
            }
            
            // Test statistics
            float min_dist = matrix.min_distance();
            float max_dist = matrix.max_distance();
            float mean_dist = matrix.mean_distance();
            
            ASSERT_GE(min_dist, 0.0f);
            ASSERT_GE(max_dist, min_dist);
            ASSERT_GE(mean_dist, min_dist);
            ASSERT_LE(mean_dist, max_dist);
        }
    }
};

} // namespace moire_test

int main() {
    auto suite = std::make_unique<moire_test::DistanceMatrixTestSuite>();
    moire_test::TestRegistry::register_suite(std::move(suite));
    moire_test::TestRegistry::run_all_tests();
    return 0;
}

