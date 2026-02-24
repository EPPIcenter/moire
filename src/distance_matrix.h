#pragma once

// Used only by cpp/tests (e.g. distance_matrix_test.cpp); not included in main package build.
// Jaccard/similarity in the main build is provided by mcmc_utils.h / UtilFunctions.

#ifndef DISTANCE_MATRIX_H_
#define DISTANCE_MATRIX_H_

#include "multivector.h"
#include <functional>
#include <string>
#include <unordered_map>
#include <memory>
#include <cmath>
#include <algorithm>
#include <numeric>

//------------------------------------------------
// Distance metric function types
template<typename T>
using DistanceFunction = std::function<T(const std::span<T const>&, const std::span<T const>&)>;

template<typename T>
using SimilarityFunction = std::function<T(const std::span<T const>&, const std::span<T const>&)>;

// Tags to distinguish between distance and similarity functions
struct DistanceTag {};
struct SimilarityTag {};

//------------------------------------------------
// Predefined distance metrics
namespace DistanceMetrics {
    
    // Jaccard distance (1 - Jaccard similarity)
    template<typename T>
    T jaccard_distance(const std::span<T const>& a, const std::span<T const>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must have the same size for Jaccard distance");
        }
        
        T intersection = 0;
        T union_sum = 0;
        
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] > 0 && b[i] > 0) {
                intersection += std::min(a[i], b[i]);
            }
            if (a[i] > 0 || b[i] > 0) {
                union_sum += std::max(a[i], b[i]);
            }
        }
        
        if (union_sum == 0) return 0;
        return 1.0 - (intersection / union_sum);
    }
    
    // Jaccard similarity
    template<typename T>
    T jaccard_similarity(const std::span<T const>& a, const std::span<T const>& b) {
        return 1.0 - jaccard_distance(a, b);
    }
    
    // Euclidean distance
    template<typename T>
    T euclidean_distance(const std::span<T const>& a, const std::span<T const>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must have the same size for Euclidean distance");
        }
        
        T sum_squares = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            T diff = a[i] - b[i];
            sum_squares += diff * diff;
        }
        
        return std::sqrt(sum_squares);
    }
    
    // Manhattan distance (L1)
    template<typename T>
    T manhattan_distance(const std::span<T const>& a, const std::span<T const>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must have the same size for Manhattan distance");
        }
        
        T sum = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            sum += std::abs(a[i] - b[i]);
        }
        
        return sum;
    }
    
    // Cosine distance (1 - cosine similarity)
    template<typename T>
    T cosine_distance(const std::span<T const>& a, const std::span<T const>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must have the same size for cosine distance");
        }
        
        T dot_product = 0;
        T norm_a = 0;
        T norm_b = 0;
        
        for (size_t i = 0; i < a.size(); ++i) {
            dot_product += a[i] * b[i];
            norm_a += a[i] * a[i];
            norm_b += b[i] * b[i];
        }
        
        norm_a = std::sqrt(norm_a);
        norm_b = std::sqrt(norm_b);
        
        if (norm_a == 0 || norm_b == 0) return 1.0;
        
        T cosine_sim = dot_product / (norm_a * norm_b);
        // Clamp to [-1, 1] to handle floating point errors
        cosine_sim = std::max(T(-1), std::min(T(1), cosine_sim));
        
        return 1.0 - cosine_sim;
    }
    
    // Cosine similarity
    template<typename T>
    T cosine_similarity(const std::span<T const>& a, const std::span<T const>& b) {
        return 1.0 - cosine_distance(a, b);
    }
    
    // Hamming distance (for binary vectors)
    template<typename T>
    T hamming_distance(const std::span<T const>& a, const std::span<T const>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must have the same size for Hamming distance");
        }
        
        T distance = 0;
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) {
                distance += 1;
            }
        }
        
        return distance;
    }
    
    // Bray-Curtis distance
    template<typename T>
    T bray_curtis_distance(const std::span<T const>& a, const std::span<T const>& b) {
        if (a.size() != b.size()) {
            throw std::invalid_argument("Vectors must have the same size for Bray-Curtis distance");
        }
        
        T sum_diff = 0;
        T sum_total = 0;
        
        for (size_t i = 0; i < a.size(); ++i) {
            sum_diff += std::abs(a[i] - b[i]);
            sum_total += a[i] + b[i];
        }
        
        if (sum_total == 0) return 0;
        return sum_diff / sum_total;
    }
}

//------------------------------------------------
// Distance matrix class
template<typename T>
class DistanceMatrix {
public:
    enum class MetricType {
        JACCARD_DISTANCE,
        JACCARD_SIMILARITY,
        EUCLIDEAN_DISTANCE,
        MANHATTAN_DISTANCE,
        COSINE_DISTANCE,
        COSINE_SIMILARITY,
        HAMMING_DISTANCE,
        BRAY_CURTIS_DISTANCE,
        CUSTOM
    };
    
    // Constructor
    DistanceMatrix(size_t num_samples, MetricType metric_type = MetricType::JACCARD_DISTANCE)
        : num_samples_(num_samples), metric_type_(metric_type) {
        
        // Initialize the distance matrix
        std::array<size_t, 2> dimensions = {num_samples, num_samples};
        distance_matrix_ = MultiVector<T, 2>(dimensions);
        
        // Set default distance function based on metric type
        set_metric(metric_type);
    }
    
    // Constructor with custom distance function
    DistanceMatrix(size_t num_samples, DistanceFunction<T> distance_func)
        : num_samples_(num_samples), metric_type_(MetricType::CUSTOM), distance_function_(distance_func) {
        
        std::array<size_t, 2> dimensions = {num_samples, num_samples};
        distance_matrix_ = MultiVector<T, 2>(dimensions);
    }
    
    // Constructor with custom similarity function
    DistanceMatrix(size_t num_samples, SimilarityFunction<T> similarity_func, bool convert_to_distance = true)
        : num_samples_(num_samples), metric_type_(MetricType::CUSTOM) {
        
        std::array<size_t, 2> dimensions = {num_samples, num_samples};
        distance_matrix_ = MultiVector<T, 2>(dimensions);
        
        if (convert_to_distance) {
            distance_function_ = [similarity_func](const std::span<T const>& a, const std::span<T const>& b) {
                return 1.0 - similarity_func(a, b);
            };
        } else {
            distance_function_ = similarity_func;
        }
    }
    
    // Constructor with custom distance function (using tag)
    template<typename F>
    DistanceMatrix(size_t num_samples, F&& distance_func, DistanceTag)
        : num_samples_(num_samples), metric_type_(MetricType::CUSTOM) {
        
        std::array<size_t, 2> dimensions = {num_samples, num_samples};
        distance_matrix_ = MultiVector<T, 2>(dimensions);
        distance_function_ = std::forward<F>(distance_func);
    }
    
    // Constructor with custom similarity function (using tag)
    template<typename F>
    DistanceMatrix(size_t num_samples, F&& similarity_func, SimilarityTag, bool convert_to_distance = true)
        : num_samples_(num_samples), metric_type_(MetricType::CUSTOM) {
        
        std::array<size_t, 2> dimensions = {num_samples, num_samples};
        distance_matrix_ = MultiVector<T, 2>(dimensions);
        
        if (convert_to_distance) {
            distance_function_ = [similarity_func = std::forward<F>(similarity_func)](const std::span<T const>& a, const std::span<T const>& b) {
                return 1.0 - similarity_func(a, b);
            };
        } else {
            distance_function_ = std::forward<F>(similarity_func);
        }
    }
    
    // Set the distance metric
    void set_metric(MetricType metric_type) {
        metric_type_ = metric_type;
        
        switch (metric_type) {
            case MetricType::JACCARD_DISTANCE:
                distance_function_ = DistanceMetrics::jaccard_distance<T>;
                break;
            case MetricType::JACCARD_SIMILARITY:
                distance_function_ = [](const std::span<T const>& a, const std::span<T const>& b) {
                    return 1.0 - DistanceMetrics::jaccard_similarity(a, b);
                };
                break;
            case MetricType::EUCLIDEAN_DISTANCE:
                distance_function_ = DistanceMetrics::euclidean_distance<T>;
                break;
            case MetricType::MANHATTAN_DISTANCE:
                distance_function_ = DistanceMetrics::manhattan_distance<T>;
                break;
            case MetricType::COSINE_DISTANCE:
                distance_function_ = DistanceMetrics::cosine_distance<T>;
                break;
            case MetricType::COSINE_SIMILARITY:
                distance_function_ = [](const std::span<T const>& a, const std::span<T const>& b) {
                    return 1.0 - DistanceMetrics::cosine_similarity(a, b);
                };
                break;
            case MetricType::HAMMING_DISTANCE:
                distance_function_ = DistanceMetrics::hamming_distance<T>;
                break;
            case MetricType::BRAY_CURTIS_DISTANCE:
                distance_function_ = DistanceMetrics::bray_curtis_distance<T>;
                break;
            case MetricType::CUSTOM:
                // Keep existing custom function
                break;
        }
    }
    
    // Set custom distance function
    void set_custom_distance(DistanceFunction<T> distance_func) {
        metric_type_ = MetricType::CUSTOM;
        distance_function_ = distance_func;
    }
    
    // Set custom similarity function
    void set_custom_similarity(SimilarityFunction<T> similarity_func, bool convert_to_distance = true) {
        metric_type_ = MetricType::CUSTOM;
        if (convert_to_distance) {
            distance_function_ = [similarity_func](const std::span<T const>& a, const std::span<T const>& b) {
                return 1.0 - similarity_func(a, b);
            };
        } else {
            distance_function_ = similarity_func;
        }
    }
    
    // Compute distance matrix from data
    template<typename DataContainer>
    void compute(const DataContainer& data) {
        if (!distance_function_) {
            throw std::runtime_error("No distance function set");
        }
        
        // Compute distances for all pairs
        for (size_t i = 0; i < num_samples_; ++i) {
            for (size_t j = 0; j < num_samples_; ++j) {
                if (i == j) {
                    distance_matrix_.at({i, j}) = 0; // Distance to self is 0
                } else if (i < j) {
                    // Compute distance only once for symmetric matrix
                    T distance = distance_function_(data[i], data[j]);
                    distance_matrix_.at({i, j}) = distance;
                    distance_matrix_.at({j, i}) = distance; // Symmetric
                }
            }
        }
    }
    
    // Compute distance matrix from data with parallel processing when enabled
    // @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    //       Falls back gracefully to sequential computation if parallel libraries unavailable.
    template<typename DataContainer>
    void compute_parallel(const DataContainer& data) {
        if (!distance_function_) {
            throw std::runtime_error("No distance function set");
        }
        
        const size_t total_elements = num_samples_ * num_samples_;
        constexpr size_t parallel_threshold = 10000; // Threshold for parallel execution
        
        // Use parallel processing for large matrices when enabled
        #ifdef MOIRE_ENABLE_PARALLEL
        if (total_elements >= parallel_threshold) {
            moire_parallel::parallel_for_2d(0, num_samples_, 0, num_samples_,
                [&](size_t i, size_t j) {
                    if (i == j) {
                        distance_matrix_.at({i, j}) = 0;
                    } else {
                        T distance = distance_function_(data[i], data[j]);
                        distance_matrix_.at({i, j}) = distance;
                    }
                });
            return;
        }
        #endif
        // Fallback to sequential computation
        compute(data);
    }
    
    // Get distance between two samples
    T get_distance(size_t i, size_t j) const {
        if (i >= num_samples_ || j >= num_samples_) {
            throw std::out_of_range("Sample indices out of range");
        }
        return distance_matrix_.at({i, j});
    }
    
    // Get similarity between two samples (1 - distance)
    T get_similarity(size_t i, size_t j) const {
        return 1.0 - get_distance(i, j);
    }
    
    // Get the entire distance matrix
    const MultiVector<T, 2>& get_matrix() const {
        return distance_matrix_;
    }
    
    // Get a row of the distance matrix
    std::vector<T> get_row(size_t i) const {
        if (i >= num_samples_) {
            throw std::out_of_range("Row index out of range");
        }
        
        std::vector<T> row(num_samples_);
        for (size_t j = 0; j < num_samples_; ++j) {
            row[j] = distance_matrix_.at({i, j});
        }
        return row;
    }
    
    // Get a column of the distance matrix
    std::vector<T> get_column(size_t j) const {
        if (j >= num_samples_) {
            throw std::out_of_range("Column index out of range");
        }
        
        std::vector<T> column(num_samples_);
        for (size_t i = 0; i < num_samples_; ++i) {
            column[i] = distance_matrix_.at({i, j});
        }
        return column;
    }
    
    // Get the number of samples
    size_t num_samples() const {
        return num_samples_;
    }
    
    // Get the metric type
    MetricType metric_type() const {
        return metric_type_;
    }
    
    // Get metric name as string
    std::string metric_name() const {
        switch (metric_type_) {
            case MetricType::JACCARD_DISTANCE: return "Jaccard Distance";
            case MetricType::JACCARD_SIMILARITY: return "Jaccard Similarity";
            case MetricType::EUCLIDEAN_DISTANCE: return "Euclidean Distance";
            case MetricType::MANHATTAN_DISTANCE: return "Manhattan Distance";
            case MetricType::COSINE_DISTANCE: return "Cosine Distance";
            case MetricType::COSINE_SIMILARITY: return "Cosine Similarity";
            case MetricType::HAMMING_DISTANCE: return "Hamming Distance";
            case MetricType::BRAY_CURTIS_DISTANCE: return "Bray-Curtis Distance";
            case MetricType::CUSTOM: return "Custom Metric";
            default: return "Unknown Metric";
        }
    }
    
    // Check if matrix is symmetric
    bool is_symmetric() const {
        for (size_t i = 0; i < num_samples_; ++i) {
            for (size_t j = i + 1; j < num_samples_; ++j) {
                if (std::abs(distance_matrix_.at({i, j}) - distance_matrix_.at({j, i})) > 1e-10) {
                    return false;
                }
            }
        }
        return true;
    }
    
    // Get minimum distance (excluding diagonal)
    T min_distance() const {
        T min_dist = std::numeric_limits<T>::max();
        for (size_t i = 0; i < num_samples_; ++i) {
            for (size_t j = 0; j < num_samples_; ++j) {
                if (i != j) {
                    min_dist = std::min(min_dist, distance_matrix_.at({i, j}));
                }
            }
        }
        return min_dist;
    }
    
    // Get maximum distance
    T max_distance() const {
        T max_dist = std::numeric_limits<T>::lowest();
        for (size_t i = 0; i < num_samples_; ++i) {
            for (size_t j = 0; j < num_samples_; ++j) {
                max_dist = std::max(max_dist, distance_matrix_.at({i, j}));
            }
        }
        return max_dist;
    }
    
    // Get mean distance (excluding diagonal)
    T mean_distance() const {
        T sum = 0;
        size_t count = 0;
        for (size_t i = 0; i < num_samples_; ++i) {
            for (size_t j = 0; j < num_samples_; ++j) {
                if (i != j) {
                    sum += distance_matrix_.at({i, j});
                    count++;
                }
            }
        }
        return count > 0 ? sum / count : 0;
    }
    
    // Clear the matrix
    void clear() {
        distance_matrix_.fill(0);
    }
    
    // Resize the matrix
    void resize(size_t new_num_samples) {
        num_samples_ = new_num_samples;
        std::array<size_t, 2> dimensions = {num_samples_, num_samples_};
        distance_matrix_.resize(dimensions);
    }

private:
    size_t num_samples_;
    MetricType metric_type_;
    DistanceFunction<T> distance_function_;
    MultiVector<T, 2> distance_matrix_;
};

#endif  // DISTANCE_MATRIX_H_
