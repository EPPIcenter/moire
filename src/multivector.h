#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <array>
#include <span>
#include <cmath>
#include <execution>
#include <ranges>

// Check if TBB is available at compile time
#if defined(__has_include) && __has_include(<tbb/parallel_for.h>) && __has_include(<tbb/blocked_range2d.h>)
#include <tbb/parallel_for.h>
#include <tbb/blocked_range2d.h>
#define TBB_AVAILABLE 1
#endif

// Check if RcppParallel is available at compile time, which includes TBB
#ifndef TBB_AVAILABLE
#if defined(__has_include) && __has_include(<RcppParallel.h>)
#include <RcppParallel.h>
#define RCPP_PARALLEL_AVAILABLE 1
// If RcppParallel is available, check if it uses TBB
#if RCPP_PARALLEL_USE_TBB
#define TBB_AVAILABLE 1
#endif
#endif
#endif

#if defined(_OPENMP)
#define OMP_AVAILABLE 1
#include <omp.h>
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define omp_get_thread_limit() 1
#define omp_get_num_procs() 1
#define omp_set_nested(a)
#define omp_set_num_threads(a)
#define omp_get_wtime() 0
#endif

#if defined(__cpp_lib_execution) && __cpp_lib_execution >= 201603
#include <execution>
#define HAS_EXECUTION 1
#endif


/// MultiVector class template
/// This class represents a multi-dimensional vector with fixed dimensions.
/// It provides methods to access elements, get sizes, and iterate over dimensions.
template <typename T, size_t N>
class MultiVector {
public:
    /// Constructor that takes dimensions
    /// @param dimensions The dimensions of the MultiVector.
    /// @throws std::invalid_argument if dimensions are empty or do not match N.
    MultiVector(const std::array<size_t, N>& dimensions) : dimensions_(dimensions) {
#ifndef NDEBUG
        if (dimensions.empty()) {
            throw std::invalid_argument("Dimensions cannot be empty");
        }

        if (N != dimensions.size()) {
            throw std::invalid_argument("Dimensions size does not match N");
        }
#endif
        // Calculate total size and strides
        size_t total_size = 1;
        for (size_t i = N; i-- > 0;) {
            strides_[i] = total_size;
            total_size *= dimensions_[i];
        }
        data_.resize(total_size, T{});
    }

    /// Constructor that takes a vector of data and dimensions
    /// @param data The data to be stored in the MultiVector.
    /// @param dimensions The dimensions of the MultiVector.
    /// @throws std::invalid_argument if dimensions are empty or do not match N.
    MultiVector(const std::vector<T>& data, const std::array<size_t, N>& dimensions) : dimensions_(dimensions) {
#ifndef NDEBUG
        if (dimensions.empty()) {
            throw std::invalid_argument("Dimensions cannot be empty");
        }
        if (N != dimensions.size()) {
            throw std::invalid_argument("Dimensions size does not match N");
        }
#endif

        size_t total_size = 1;
        for (size_t i = N; i-- > 0;) {
            strides_[i] = total_size;
            total_size *= dimensions_[i];
        }
        if (data.size() != total_size) {
            throw std::invalid_argument("Data size does not match total size");
        }
        data_ = data;
    }

    MultiVector() {
        dimensions_ = {0};
        data_ = {};
        strides_ = {};

    }

    void resize(const std::array<size_t, N>& dimensions) {
#ifndef NDEBUG
        if (dimensions.empty()) {
            throw std::invalid_argument("Dimensions cannot be empty");
        }
        if (N != dimensions.size()) {
            throw std::invalid_argument("Dimensions size does not match N");
        }
#endif
        dimensions_ = dimensions;
        size_t total_size = 1;
        for (size_t i = N; i-- > 0;) {
            strides_[i] = total_size;
            total_size *= dimensions_[i];
        }
        data_.resize(total_size, T{});
    }
    /// Fill the innermost dimension with a single value
    /// @param indices The indices of the outer dimensions.
    /// @param value The value to fill the innermost dimension with.
    void inner_fill(const std::array<size_t, N - 1>& indices, const T& value) {
        auto [begin, end] = inner_iterators(indices);
        std::fill(begin, end, value);
    }

    /// Fill the innermost dimension with a range of values
    /// @param indices The indices of the outer dimensions.
    /// @param values The values to fill the innermost dimension with.
    void inner_fill(const std::array<size_t, N - 1>& indices, const std::span<T const> values) {
#ifndef NDEBUG
        if (values.size() > dimensions_.at(indices.back())) {
            throw std::invalid_argument("(" + std::to_string(N) + "D) Incorrect number of values, expected " + std::to_string(dimensions_.at(indices.back())) + " but got " + std::to_string(values.size()));
        }
#endif
        auto [begin, end] = inner_iterators(indices);
        std::copy(values.begin(), values.end(), begin);
    }

    /// Access element at a given multi-dimensional index
    /// @param indices The indices of the element to access.
    /// @return A reference to the element at the specified indices.
    /// @throws std::invalid_argument if the number of indices is incorrect.
    T& at(const std::array<size_t, N>& indices) {
#ifndef NDEBUG
        if (indices.size() != dimensions_.size()) {
            throw std::invalid_argument("Incorrect number of indices");
        }
#endif
        return data_.at(calculate_index(indices));
    }

    const T& at(const std::array<size_t, N>& indices) const {
#ifndef NDEBUG
        if (indices.size() != dimensions_.size()) {
            throw std::invalid_argument("Incorrect number of indices");
        }
#endif
        return data_.at(calculate_index(indices));
    }

    /// Get the size of a specific dimension
    /// @param dimension The dimension to query.
    /// @return The size of the specified dimension.
    /// @throws std::out_of_range if the dimension is out of range.
    size_t size(size_t dimension) const {
#ifndef NDEBUG
        if (dimension >= dimensions_.size()) {
            throw std::out_of_range("Dimension out of range");
        }
#endif
        return dimensions_[dimension];
    }

    /// Get the total number of elements
    /// @return The total number of elements in the MultiVector.
    size_t total_size() const {
        return data_.size();
    }

    /// Iterator for the innermost dimension
    /// @param outer_indices The indices of the outer dimensions.
    /// @return An iterator to the beginning of the innermost dimension.
    /// @throws std::invalid_argument if the number of outer indices is incorrect.
    typename std::vector<T>::iterator inner_begin(const std::array<size_t, N-1>& outer_indices) {
#ifndef NDEBUG
        if (outer_indices.size() != dimensions_.size() - 1) {
            throw std::invalid_argument("Incorrect number of outer indices");
        }
#endif
        std::array<size_t, N> full_indices{0};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        size_t start_index = calculate_index(full_indices);
        return data_.begin() + start_index;
    }

    typename std::vector<T>::iterator inner_end(const std::array<size_t, N-1>& outer_indices) {
#ifndef NDEBUG
        if (outer_indices.size() != dimensions_.size() - 1) {
            throw std::invalid_argument("Incorrect number of outer indices");
        }
#endif
        std::array<size_t, N> full_indices{0};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        size_t start_index = calculate_index(full_indices);
        return data_.begin() + start_index + dimensions_.back();
    }

    /// Const iterator for the innermost dimension
    /// @param outer_indices The indices of the outer dimensions.
    /// @return A const iterator to the beginning of the innermost dimension.
    /// @throws std::invalid_argument if the number of outer indices is incorrect.
    typename std::vector<T>::const_iterator inner_begin(const std::array<size_t, N-1>& outer_indices) const {
#ifndef NDEBUG
        if (outer_indices.size() != dimensions_.size() - 1) {
            throw std::invalid_argument("Incorrect number of outer indices");
        }
#endif
        std::array<size_t, N> full_indices{0};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        size_t start_index = calculate_index(full_indices);
        return data_.cbegin() + start_index;
    }

    typename std::vector<T>::const_iterator inner_end(const std::array<size_t, N-1>& outer_indices) const {
#ifndef NDEBUG
        if (outer_indices.size() != dimensions_.size() - 1) {
            throw std::invalid_argument("Incorrect number of outer indices");
        }
#endif
        std::array<size_t, N> full_indices{0};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        size_t start_index = calculate_index(full_indices);
        return data_.cbegin() + start_index + dimensions_.back();
    }

    std::pair<typename std::vector<T>::const_iterator, typename std::vector<T>::const_iterator> inner_iterators(const std::array<size_t, N - 1>& outer_indices) const {
        auto [start_index, end_index] = calculate_start_end_indices(outer_indices);
        return {data_.cbegin() + start_index, data_.cbegin() + end_index};
    }

    std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator> inner_iterators(const std::array<size_t, N - 1>& outer_indices) {
        auto [start_index, end_index] = calculate_start_end_indices(outer_indices);
        return {data_.begin() + start_index, data_.begin() + end_index};
    }

    /// Clear the MultiVector
    /// @note This will clear the data contained in the MultiVector.
    void clear() {
        data_.clear();
    }

    std::array<size_t, N> dimensions() const {
        return dimensions_;
    }

    const std::array<size_t, N>& strides() const {
        return strides_;
    }

    const std::vector<T>& data() const {
        return data_;
    }

    /// Reduce the innermost dimension using a binary operation, returning a MultiVector of dimension N-1
    /// @param binary_op The binary operation to use for reduction
    /// @param init The initial value for the reduction
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A MultiVector<T, N-1> containing the reduced values
    template<typename BinaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    auto reduce(BinaryOp binary_op, const auto& init, ExecutionPolicy policy = std::execution::seq) const {
        using ResultType = decltype(binary_op(init, std::declval<T>()));
        // Create output MultiVector with dimensions excluding the last dimension
        std::array<size_t, N-1> reduced_dims;
        std::copy(dimensions_.begin(), dimensions_.end() - 1, reduced_dims.begin());
        MultiVector<ResultType, N-1> result(reduced_dims);
        
        // Optimized implementation for 3D case (most common in benchmarks)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
            // Use a more direct approach for 3D case
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    
                    // Perform reduction
                    result.at({i, j}) = std::reduce(policy, begin, end, init, binary_op);
                }
            }
        }
        // Optimized implementation for 2D case
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
            // Use a direct approach for 2D case
            for (size_t i = 0; i < dim0; ++i) {
                // Calculate the start index for this row
                const size_t start_idx = i * stride0;
                
                // Get iterators for the innermost dimension
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                
                // Perform reduction
                result.at({i}) = std::reduce(policy, begin, end, init, binary_op);
            }
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto reduce_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Perform reduction for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    result.at(indices) = std::reduce(policy, begin, end, init, binary_op);
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            reduce_recursive(reduce_recursive, indices, 0);
        }
        
        return result;
    }

    /// Optimized parallel reduce using TBB for better performance if available, otherwise falls back to standard library
    /// @param binary_op The binary operation to use for reduction
    /// @param init The initial value for the reduction
    /// @return A MultiVector<T, N-1> containing the reduced values
    template<typename BinaryOp>
#if !defined(TBB_AVAILABLE) && !defined(HAS_EXECUTION) && !defined(OMP_AVAILABLE)
    [[deprecated("Consider enabling TBB, C++17 execution policies, or OpenMP for better performance.")]]
#endif
    auto parallel_reduce(BinaryOp binary_op, const auto& init) const {
        using ResultType = decltype(binary_op(init, std::declval<T>()));
        // Create output MultiVector with dimensions excluding the last dimension
        std::array<size_t, N-1> reduced_dims;
        std::copy(dimensions_.begin(), dimensions_.end() - 1, reduced_dims.begin());
        MultiVector<ResultType, N-1> result(reduced_dims);
        
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
#if defined(TBB_AVAILABLE)
            // Use TBB's parallel_for with a 2D blocked range for better load balancing
            tbb::parallel_for(
                tbb::blocked_range2d<size_t>(0, dim0, 0, dim1),
                [&](const tbb::blocked_range2d<size_t>& range) {
                    for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
                        for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
                            // Calculate the start index for this 2D slice
                            const size_t start_idx = i * stride0 + j * stride1;
                            
                            // Get iterators for the innermost dimension
                            auto begin = data_.begin() + start_idx;
                            auto end = begin + dim2;
                            
                            // Perform reduction
                            result.at({i, j}) = std::reduce(std::execution::seq, begin, end, init, binary_op);
                        }
                    }
                }
            );
#elif defined(HAS_EXECUTION)
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0 * dim1);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::par_unseq, 
                indices.begin(), 
                indices.end(),
                [&](size_t idx) {
                    const size_t i = idx / dim1;
                    const size_t j = idx % dim1;
                    
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    
                    // Perform reduction
                    result.at({i, j}) = std::reduce(std::execution::seq, begin, end, init, binary_op);
                }
            );
#elif defined(OMP_AVAILABLE)
            #pragma omp parallel for collapse(2)
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    
                    // Perform reduction
                    result.at({i, j}) = std::reduce(std::execution::seq, begin, end, init, binary_op);
                }
            }
#else
            // Sequential fallback
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    
                    // Perform reduction
                    result.at({i, j}) = std::reduce(std::execution::seq, begin, end, init, binary_op);
                }
            }
#endif
        }
        // Specialized implementation for 2D case
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
#if defined(TBB_AVAILABLE)
            // Use TBB's parallel_for for 1D case
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, dim0),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i != range.end(); ++i) {
                        // Calculate the start index for this row
                        const size_t start_idx = i * stride0;
                        
                        // Get iterators for the innermost dimension
                        auto begin = data_.begin() + start_idx;
                        auto end = begin + dim1;
                        
                        // Perform reduction
                        result.at({i}) = std::reduce(std::execution::seq, begin, end, init, binary_op);
                    }
                }
            );
#elif defined(HAS_EXECUTION)
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::par_unseq, 
                indices.begin(), 
                indices.end(),
                [&](size_t i) {
                    // Calculate the start index for this row
                    const size_t start_idx = i * stride0;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim1;
                    
                    // Perform reduction
                    result.at({i}) = std::reduce(std::execution::seq, begin, end, init, binary_op);
                }
            );
#elif defined(OMP_AVAILABLE)
            #pragma omp parallel for
            for (size_t i = 0; i < dim0; ++i) {
                // Calculate the start index for this row
                const size_t start_idx = i * stride0;
                
                // Get iterators for the innermost dimension
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                
                // Perform reduction
                result.at({i}) = std::reduce(std::execution::seq, begin, end, init, binary_op);
            }
#else
            // Sequential fallback
            for (size_t i = 0; i < dim0; ++i) {
                // Calculate the start index for this row
                const size_t start_idx = i * stride0;
                
                // Get iterators for the innermost dimension
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                
                // Perform reduction
                result.at({i}) = std::reduce(std::execution::seq, begin, end, init, binary_op);
            }
#endif
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto reduce_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Perform reduction for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    result.at(indices) = std::reduce(std::execution::seq, begin, end, init, binary_op);
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            reduce_recursive(reduce_recursive, indices, 0);
        }
        
        return result;
    }

    
    /**
     * @brief Fully reduces the multivector down to a scalar value using a binary operation
     * 
     * This function applies a binary operation to all elements in the multivector,
     * reducing them to a single scalar value.
     * 
     * @tparam BinaryOp The type of the binary operation to apply
     * @tparam ExecutionPolicy The execution policy for parallelization
     * @param binary_op The binary operation to apply (e.g., std::plus, std::multiplies)
     * @param policy The execution policy (default: sequential)
     * @return A scalar value containing the fully reduced result
     */
    template<typename BinaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    T full_reduce(BinaryOp binary_op, ExecutionPolicy policy = std::execution::seq) const {
        return std::reduce(policy, data_.begin(), data_.end(), T{}, binary_op);
    }

    /**
     * @brief Fully reduces the multivector down to a scalar value in parallel
     * 
     * This function applies a binary operation to all elements in the multivector,
     * reducing them to a single scalar value using parallel processing when TBB is available.
     * 
     * @tparam BinaryOp The type of the binary operation to apply
     * @param binary_op The binary operation to apply (e.g., std::plus, std::multiplies)
     * @return A scalar value containing the fully reduced result
     */
    template<typename BinaryOp>
    T parallel_full_reduce(BinaryOp binary_op) const {
        #if TBB_AVAILABLE
        return tbb::parallel_reduce(
            tbb::blocked_range<typename std::vector<T>::const_iterator>(data_.begin(), data_.end()),
            T{},
            [&binary_op](const tbb::blocked_range<typename std::vector<T>::const_iterator>& range, T init) {
                return std::reduce(std::execution::seq, range.begin(), range.end(), init, binary_op);
            }
        );
        #else
        return std::reduce(std::execution::seq, data_.begin(), data_.end(), init, binary_op);
        #endif
    }

    /**
     * @brief Computes the sum of all elements in the multivector
     * 
     * @tparam ExecutionPolicy The execution policy for parallelization
     * @param policy The execution policy (default: sequential)
     * @return The sum of all elements
     */
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T full_sum(ExecutionPolicy policy = std::execution::seq) const {
        return full_reduce(std::plus<T>(), policy);
    }

    /**
     * @brief Computes the sum of all elements in the multivector in parallel
     * 
     * @return The sum of all elements
     */
    T parallel_full_sum() const {
        return parallel_full_reduce(std::plus<T>());
    }

    /**
     * @brief Computes the product of all elements in the multivector
     * 
     * @tparam ExecutionPolicy The execution policy for parallelization
     * @param policy The execution policy (default: sequential)
     * @return The product of all elements
     */
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T full_product(ExecutionPolicy policy = std::execution::seq) const {
        return full_reduce(std::multiplies<T>(), policy);
    }

    /**
     * @brief Computes the product of all elements in the multivector in parallel
     * 
     * @return The product of all elements
     */
    T parallel_full_product() const {
        return parallel_full_reduce(std::multiplies<T>());
    }

    /**
     * @brief Finds the minimum element in the multivector
     * 
     * @tparam ExecutionPolicy The execution policy for parallelization
     * @param policy The execution policy (default: sequential)
     * @return The minimum element
     */
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T full_min(ExecutionPolicy policy = std::execution::seq) const {
        return full_reduce(std::min<T>(), policy);
    }

    /**
     * @brief Finds the minimum element in the multivector in parallel
     * 
     * @return The minimum element
     */
    T parallel_full_min() const {
        return parallel_full_reduce(std::min<T>());
    }

    /**
     * @brief Finds the maximum element in the multivector
     * 
     * @tparam ExecutionPolicy The execution policy for parallelization
     * @param policy The execution policy (default: sequential)
     * @return The maximum element
     */
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T full_max(ExecutionPolicy policy = std::execution::seq) const {
        return full_reduce(std::max<T>(), policy);
    }

    /**
     * @brief Finds the maximum element in the multivector in parallel
     * 
     * @return The maximum element
     */
    T parallel_full_max() const {
        return parallel_full_reduce(std::max<T>());
    }

    // Convenience methods for common reductions with optional execution policy
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N-1> sum(ExecutionPolicy policy = std::execution::seq) const {
        return reduce(std::plus<T>(), T{}, policy);
    }

    // Optimized parallel sum method
    MultiVector<T, N-1> parallel_sum() const {
        return parallel_reduce(std::plus<T>(), T{});
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N-1> product(ExecutionPolicy policy = std::execution::seq) const {
        return reduce(std::multiplies<T>(), T{1}, policy);
    }

    // Optimized parallel product method
    MultiVector<T, N-1> parallel_product() const {
        return parallel_reduce(std::multiplies<T>(), T{1});
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N-1> max(ExecutionPolicy policy = std::execution::seq) const {
        return reduce([](const T& a, const T& b) { return std::max(a, b); }, std::numeric_limits<T>::lowest(), policy);
    }

    // Optimized parallel max method
    MultiVector<T, N-1> parallel_max() const {
        return parallel_reduce([](const T& a, const T& b) { return std::max(a, b); },std::numeric_limits<T>::lowest());
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N-1> min(ExecutionPolicy policy = std::execution::seq) const {
        return reduce([](const T& a, const T& b) { return std::min(a, b); }, std::numeric_limits<T>::max(), policy);
    }

    // Optimized parallel min method
    MultiVector<T, N-1> parallel_min() const {
        return parallel_reduce([](const T& a, const T& b) { return std::min(a, b); },
                              std::numeric_limits<T>::max());
    }

    /// Compute the log-sum-exp reduction over the innermost dimension
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A MultiVector<T, N-1> containing the log-sum-exp values
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N-1> logsumexp(ExecutionPolicy policy = std::execution::seq) const {
        // Create output MultiVector with dimensions excluding the last dimension
        std::array<size_t, N-1> reduced_dims;
        std::copy(dimensions_.begin(), dimensions_.end() - 1, reduced_dims.begin());
        MultiVector<T, N-1> result(reduced_dims);
        
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
            // Use a more direct approach for 3D case
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    
                    // Find the maximum value for numerical stability
                    T max_val = *std::max_element(policy, begin, end);
                    
                    // Compute exp(x - max) and sum
                    T sum = T{};
                    for (auto it = begin; it != end; ++it) {
                        sum += std::exp(*it - max_val);
                    }
                    
                    // Compute log(sum) + max
                    result.at({i, j}) = std::log(sum) + max_val;
                }
            }
        }
        // Optimized implementation for 2D case
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
            // Use a direct approach for 2D case
            for (size_t i = 0; i < dim0; ++i) {
                // Calculate the start index for this row
                const size_t start_idx = i * stride0;
                
                // Get iterators for the innermost dimension
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                
                // Find the maximum value for numerical stability
                T max_val = *std::max_element(policy, begin, end);
                
                // Compute exp(x - max) and sum
                T sum = T{};
                for (auto it = begin; it != end; ++it) {
                    sum += std::exp(*it - max_val);
                }
                
                // Compute log(sum) + max
                result.at({i}) = std::log(sum) + max_val;
            }
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto logsumexp_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Perform logsumexp for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    
                    // Find the maximum value for numerical stability
                    T max_val = *std::max_element(policy, begin, end);
                    
                    // Compute exp(x - max) and sum
                    T sum = T{};
                    for (auto it = begin; it != end; ++it) {
                        sum += std::exp(*it - max_val);
                    }
                    
                    // Compute log(sum) + max
                    result.at(indices) = std::log(sum) + max_val;
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            logsumexp_recursive(logsumexp_recursive, indices, 0);
        }
        
        return result;
    }

    /// Compute the log-sum-exp reduction over the innermost dimension in parallel
    /// @return A MultiVector<T, N-1> containing the log-sum-exp values
    MultiVector<T, N-1> parallel_logsumexp() const {
        // Create output MultiVector with dimensions excluding the last dimension
        std::array<size_t, N-1> reduced_dims;
        std::copy(dimensions_.begin(), dimensions_.end() - 1, reduced_dims.begin());
        MultiVector<T, N-1> result(reduced_dims);
        
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for with a 2D blocked range for better load balancing
            tbb::parallel_for(
                tbb::blocked_range2d<size_t>(0, dim0, 0, dim1),
                [&](const tbb::blocked_range2d<size_t>& range) {
                    for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
                        for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
                            // Calculate the start index for this 2D slice
                            const size_t start_idx = i * stride0 + j * stride1;
                            
                            // Get iterators for the innermost dimension
                            auto begin = data_.begin() + start_idx;
                            auto end = begin + dim2;
                            
                            // Find the maximum value for numerical stability
                            T max_val = *std::max_element(std::execution::seq, begin, end);
                            
                            // Compute exp(x - max) and sum
                            T sum = T{};
                            for (auto it = begin; it != end; ++it) {
                                sum += std::exp(*it - max_val);
                            }
                            
                            // Compute log(sum) + max
                            result.at({i, j}) = std::log(sum) + max_val;
                        }
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0 * dim1);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::seq, 
                indices.begin(), 
                indices.end(),
                [&](size_t idx) {
                    const size_t i = idx / dim1;
                    const size_t j = idx % dim1;
                    
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    
                    // Find the maximum value for numerical stability
                    T max_val = *std::max_element(std::execution::seq, begin, end);
                    
                    // Compute exp(x - max) and sum
                    T sum = T{};
                    for (auto it = begin; it != end; ++it) {
                        sum += std::exp(*it - max_val);
                    }
                    
                    // Compute log(sum) + max
                    result.at({i, j}) = std::log(sum) + max_val;
                }
            );
#endif
        }
        // Specialized implementation for 2D case
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for for 1D case
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, dim0),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i != range.end(); ++i) {
                        // Calculate the start index for this row
                        const size_t start_idx = i * stride0;
                        
                        // Get iterators for the innermost dimension
                        auto begin = data_.begin() + start_idx;
                        auto end = begin + dim1;
                        
                        // Find the maximum value for numerical stability
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        
                        // Compute exp(x - max) and sum
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        
                        // Compute log(sum) + max
                        result.at({i}) = std::log(sum) + max_val;
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::seq, 
                indices.begin(), 
                indices.end(),
                [&](size_t i) {
                    // Calculate the start index for this row
                    const size_t start_idx = i * stride0;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim1;
                    
                    // Find the maximum value for numerical stability
                    T max_val = *std::max_element(std::execution::seq, begin, end);
                    
                    // Compute exp(x - max) and sum
                    T sum = T{};
                    for (auto it = begin; it != end; ++it) {
                        sum += std::exp(*it - max_val);
                    }
                    
                    // Compute log(sum) + max
                    result.at({i}) = std::log(sum) + max_val;
                }
            );
#endif
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto logsumexp_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Perform logsumexp for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    
                    // Find the maximum value for numerical stability
                    T max_val = *std::max_element(std::execution::seq, begin, end);
                    
                    // Compute exp(x - max) and sum
                    T sum = T{};
                    for (auto it = begin; it != end; ++it) {
                        sum += std::exp(*it - max_val);
                    }
                    
                    // Compute log(sum) + max
                    result.at(indices) = std::log(sum) + max_val;
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            logsumexp_recursive(logsumexp_recursive, indices, 0);
        }
        
        return result;
    }

    /// Apply a transformation to each inner slice of the MultiVector
    /// @param op The binary operation to apply to each element in the inner slice
    /// @param values A span of values to apply to each inner slice
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A new MultiVector with the transformed values
    template<typename BinaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> transform(BinaryOp op, const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        MultiVector<T, N> result(dimensions_);
        // Calculate the number of inner slices
        size_t num_inner_slices = 1;
        for (size_t i = 0; i < N - 1; ++i) {
            num_inner_slices *= dimensions_[i];
        }
        
#ifndef NDEBUG
        if (values.size() != num_inner_slices) {
            throw std::invalid_argument("Values span size does not match the number of inner slices");
        }
#endif
        
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    const size_t slice_idx = i * dim1 + j;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    // Apply transformation to each element in the inner slice
                    std::transform(policy, begin, end, result_begin, 
                                  [&](const T& x) { return op(x, values[slice_idx]); });
                }
            }
        }
        // Optimized implementation for 2D case
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
            for (size_t i = 0; i < dim0; ++i) {
                // Calculate the start index for this row
                const size_t start_idx = i * stride0;
                
                // Get iterators for the innermost dimension
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                auto result_begin = result.data_.begin() + start_idx;
                
                // Apply transformation to each element in the inner slice
                std::transform(policy, begin, end, result_begin, 
                              [&](const T& x) { return op(x, values[i]); });
            }
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim, size_t& slice_idx) -> void {
                if (dim == N-1) {
                    // Apply transformation for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    std::transform(policy, begin, end, result_begin, 
                                  [&](const T& x) { return op(x, values[slice_idx]); });
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1, slice_idx);
                }
            };

            std::array<size_t, N-1> indices{};
            size_t slice_idx = 0;
            transform_recursive(transform_recursive, indices, 0, slice_idx);
        }
        return result;
    }

    /// Apply a transformation to each inner slice of the MultiVector in parallel
    /// @param op The binary operation to apply to each element in the inner slice
    /// @param values A span of values to apply to each inner slice
    /// @return A new MultiVector with the transformed values
    template<typename BinaryOp>
    MultiVector<T, N> parallel_transform(BinaryOp op, const std::span<T const> values) const {
        MultiVector<T, N> result(dimensions_);
        // Calculate the number of inner slices
        size_t num_inner_slices = 1;
        for (size_t i = 0; i < N - 1; ++i) {
            num_inner_slices *= dimensions_[i];
        }
        
#ifndef NDEBUG
        if (values.size() != num_inner_slices) {
            throw std::invalid_argument("Values span size does not match the number of inner slices");
        }
#endif
        
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for with a 2D blocked range for better load balancing
            tbb::parallel_for(
                tbb::blocked_range2d<size_t>(0, dim0, 0, dim1),
                [&](const tbb::blocked_range2d<size_t>& range) {
                    for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
                        for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
                            // Calculate the start index for this 2D slice
                            const size_t start_idx = i * stride0 + j * stride1;
                            const size_t slice_idx = i * dim1 + j;
                            
                            // Get iterators for the innermost dimension
                            auto begin = data_.begin() + start_idx;
                            auto end = begin + dim2;
                            auto result_begin = result.data_.begin() + start_idx;
                            
                            // Apply transformation to each element in the inner slice
                            std::transform(std::execution::seq, begin, end, result_begin, 
                                          [&](const T& x) { return op(x, values[slice_idx]); });
                        }
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0 * dim1);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::seq, 
                indices.begin(), 
                indices.end(),
                [&](size_t idx) {
                    const size_t i = idx / dim1;
                    const size_t j = idx % dim1;
                    
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    const size_t slice_idx = i * dim1 + j;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    // Apply transformation to each element in the inner slice
                    std::transform(std::execution::seq, begin, end, result_begin, 
                                  [&](const T& x) { return op(x, values[slice_idx]); });
                }
            );
#endif
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for for 1D case
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, dim0),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i != range.end(); ++i) {
                        // Calculate the start index for this row
                        const size_t start_idx = i * stride0;
                        
                        // Get iterators for the innermost dimension
                        auto begin = data_.begin() + start_idx;
                        auto end = begin + dim1;
                        auto result_begin = result.data_.begin() + start_idx;
                        
                        // Apply transformation to each element in the inner slice
                        std::transform(std::execution::seq, begin, end, result_begin, 
                                      [&](const T& x) { return op(x, values[i]); });
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::seq, 
                indices.begin(), 
                indices.end(),
                [&](size_t i) {
                    // Calculate the start index for this row
                    const size_t start_idx = i * stride0;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim1;
                    
                    // Apply transformation to each element in the inner slice
                    std::transform(std::execution::seq, begin, end, begin, 
                                  [&](const T& x) { return op(x, values[i]); });
                }
            );
#endif
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim, size_t& slice_idx) -> void {
                if (dim == N-1) {
                    // Apply transformation for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    std::transform(std::execution::seq, begin, end, begin, 
                                  [&](const T& x) { return op(x, values[slice_idx]); });
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1, slice_idx);
                }
            };

            std::array<size_t, N-1> indices{};
            size_t slice_idx = 0;
            transform_recursive(transform_recursive, indices, 0, slice_idx);
        }
    }

    

    /// Add a value to each inner slice of the MultiVector
    /// @param values A span of values to add to each inner slice
    /// @param execution_policy Optional execution policy for parallel execution
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> add(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::plus<T>(), values, policy);
    }

    /// Add a value to each inner slice of the MultiVector in parallel
    /// @param values A span of values to add to each inner slice
    MultiVector<T, N> parallel_add(const std::span<T const> values) const {
        return parallel_transform(std::plus<T>(), values);
    }

    /// Subtract a value from each inner slice of the MultiVector
    /// @param values A span of values to subtract from each inner slice
    /// @param execution_policy Optional execution policy for parallel execution
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> subtract(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::minus<T>(), values, policy);
    }

    /// Subtract a value from each inner slice of the MultiVector in parallel
    /// @param values A span of values to subtract from each inner slice
    MultiVector<T, N> parallel_subtract(const std::span<T const> values) const {
        return parallel_transform(std::minus<T>(), values);
    }

    /// Multiply each inner slice of the MultiVector by a value
    /// @param values A span of values to multiply each inner slice by
    /// @param execution_policy Optional execution policy for parallel execution
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> multiply(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::multiplies<T>(), values, policy);
    }

    /// Multiply each inner slice of the MultiVector by a value in parallel
    /// @param values A span of values to multiply each inner slice by
    MultiVector<T, N> parallel_multiply(const std::span<T const> values) const {
        return parallel_transform(std::multiplies<T>(), values);
    }

    /// Divide each inner slice of the MultiVector by a value
    /// @param values A span of values to divide each inner slice by
    /// @param execution_policy Optional execution policy for parallel execution
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> divide(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::divides<T>(), values, policy);
    }

    /// Divide each inner slice of the MultiVector by a value in parallel
    /// @param values A span of values to divide each inner slice by
    MultiVector<T, N> parallel_divide(const std::span<T const> values) const {
        return parallel_transform(std::divides<T>(), values);
    }

    /// Apply a transformation to each element within each inner slice of the MultiVector
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element position across all inner slices
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A new MultiVector with the transformed values
    template<typename BinaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> element_transform(BinaryOp op, const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
#ifndef NDEBUG
        if (values.size() != dimensions_.back()) {
            throw std::invalid_argument("Values span size does not match the inner dimension size");
        }
#endif
        
        MultiVector<T, N> result(dimensions_);
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Apply transformation to each element position
                    for (size_t k = 0; k < dim2; ++k) {
                        result.data_[start_idx + k] = op(data_[start_idx + k], values[k]);
                    }
                }
            }
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
            for (size_t i = 0; i < dim0; ++i) {
                // Calculate the start index for this row
                const size_t start_idx = i * stride0;
                
                // Apply transformation to each element position
                for (size_t j = 0; j < dim1; ++j) {
                    result.data_[start_idx + j] = op(data_[start_idx + j], values[j]);
                }
            }
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Apply transformation for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    size_t pos = 0;
                    for (auto it = begin; it != end; ++it, ++pos) {
                        *(result_begin + pos) = op(*it, values[pos]);
                    }
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            transform_recursive(transform_recursive, indices, 0);
        }
        return result;
    }

    /// Apply a transformation to each element within each inner slice of the MultiVector in parallel
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element position across all inner slices
    /// @return A new MultiVector with the transformed values
    template<typename BinaryOp>
    MultiVector<T, N> parallel_element_transform(BinaryOp op, const std::span<T const> values) const {
#ifndef NDEBUG
        if (values.size() != dimensions_.back()) {
            throw std::invalid_argument("Values span size does not match the inner dimension size");
        }
#endif
        
        MultiVector<T, N> result(dimensions_);
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for with a 2D blocked range for better load balancing
            tbb::parallel_for(
                tbb::blocked_range2d<size_t>(0, dim0, 0, dim1),
                [&](const tbb::blocked_range2d<size_t>& range) {
                    for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
                        for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
                            // Calculate the start index for this 2D slice
                            const size_t start_idx = i * stride0 + j * stride1;
                            
                            // Apply transformation to each element position
                            for (size_t k = 0; k < dim2; ++k) {
                                result.data_[start_idx + k] = op(data_[start_idx + k], values[k]);
                            }
                        }
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0 * dim1);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::seq, 
                indices.begin(), 
                indices.end(),
                [&](size_t idx) {
                    const size_t i = idx / dim1;
                    const size_t j = idx % dim1;
                    
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Apply transformation to each element position
                    for (size_t k = 0; k < dim2; ++k) {
                        result.data_[start_idx + k] = op(data_[start_idx + k], values[k]);
                    }
                }
            );
#endif
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for for 1D case
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, dim0),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i != range.end(); ++i) {
                        // Calculate the start index for this row
                        const size_t start_idx = i * stride0;
                        
                        // Apply transformation to each element position
                        for (size_t j = 0; j < dim1; ++j) {
                            result.data_[start_idx + j] = op(data_[start_idx + j], values[j]);
                        }
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::seq, 
                indices.begin(), 
                indices.end(),
                [&](size_t i) {
                    // Calculate the start index for this row
                    const size_t start_idx = i * stride0;
                    
                    // Apply transformation to each element position
                    for (size_t j = 0; j < dim1; ++j) {
                        result.data_[start_idx + j] = op(data_[start_idx + j], values[j]);
                    }
                }
            );
#endif
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Apply transformation for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    size_t pos = 0;
                    for (auto it = begin; it != end; ++it, ++pos) {
                        *(result_begin + pos) = op(*it, values[pos]);
                    }
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            transform_recursive(transform_recursive, indices, 0);
        }
        return result;
    }

    // Convenience methods for element-level operations
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> element_add(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return element_transform(std::plus<T>(), values, policy);
    }

    MultiVector<T, N> parallel_element_add(const std::span<T const> values) const {
        return parallel_element_transform(std::plus<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> element_subtract(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return element_transform(std::minus<T>(), values, policy);
    }

    MultiVector<T, N> parallel_element_subtract(const std::span<T const> values) const {
        return parallel_element_transform(std::minus<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> element_multiply(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return element_transform(std::multiplies<T>(), values, policy);
    }

    MultiVector<T, N> parallel_element_multiply(const std::span<T const> values) const {
        return parallel_element_transform(std::multiplies<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> element_divide(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return element_transform(std::divides<T>(), values, policy);
    }

    MultiVector<T, N> parallel_element_divide(const std::span<T const> values) const {
        return parallel_element_transform(std::divides<T>(), values);
    }

    /// Check if two multivectors have the same dimensions
    /// @param other The other multivector to compare with
    /// @return true if dimensions match, false otherwise
    bool has_same_dimensions(const MultiVector<T, N>& other) const {
        return dimensions_ == other.dimensions_;
    }

    /// Elementwise addition of two multivectors
    /// @param other The other multivector to add
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
    MultiVector<T, N> operator+(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
        
#if TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = data_[i] + other.data_[i];
                }
            });
#elif HAS_EXECUTION
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::plus<T>());
#else
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] + other.data_[i];
        }
#endif
        return result;
    }

    /// Elementwise subtraction of two multivectors
    /// @param other The other multivector to subtract
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
    MultiVector<T, N> operator-(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
        
#if TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = data_[i] - other.data_[i];
                }
            });
#elif HAS_EXECUTION
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::minus<T>());
#else
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] - other.data_[i];
        }
#endif
        return result;
    }

    /// Elementwise multiplication of two multivectors
    /// @param other The other multivector to multiply with
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
    MultiVector<T, N> operator*(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
        
#if TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = data_[i] * other.data_[i];
                }
            });
#elif HAS_EXECUTION
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::multiplies<T>());
#else
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] * other.data_[i];
        }
#endif
        return result;
    }

    /// Elementwise division of two multivectors
    /// @param other The other multivector to divide by
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
    MultiVector<T, N> operator/(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
        
#if TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = data_[i] / other.data_[i];
                }
            });
#elif HAS_EXECUTION
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::divides<T>());
#else
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] / other.data_[i];
        }
#endif
        return result;
    }

    /// In-place elementwise addition
    /// @param other The other multivector to add
    /// @return Reference to this multivector
    /// @throws std::invalid_argument if dimensions don't match
    MultiVector<T, N>& operator+=(const MultiVector<T, N>& other) {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
#if TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    data_[i] += other.data_[i];
                }
            });
#elif HAS_EXECUTION
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      data_.begin(), std::plus<T>());
#else
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] += other.data_[i];
        }
#endif
        return *this;
    }

    /// In-place elementwise subtraction
    /// @param other The other multivector to subtract
    /// @return Reference to this multivector
    /// @throws std::invalid_argument if dimensions don't match
    MultiVector<T, N>& operator-=(const MultiVector<T, N>& other) {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
#if TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    data_[i] -= other.data_[i];
                }
            });
#elif HAS_EXECUTION
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      data_.begin(), std::minus<T>());
#else
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] -= other.data_[i];
        }
#endif
        return *this;
    }

    /// In-place elementwise multiplication
    /// @param other The other multivector to multiply with
    /// @return Reference to this multivector
    /// @throws std::invalid_argument if dimensions don't match
    MultiVector<T, N>& operator*=(const MultiVector<T, N>& other) {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
#if TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    data_[i] *= other.data_[i];
                }
            });
#elif HAS_EXECUTION
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      data_.begin(), std::multiplies<T>());
#else
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] *= other.data_[i];
        }
#endif
        return *this;
    }

    /// In-place elementwise division
    /// @param other The other multivector to divide by
    /// @return Reference to this multivector
    /// @throws std::invalid_argument if dimensions don't match
    MultiVector<T, N>& operator/=(const MultiVector<T, N>& other) {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
#if TBB_AVAILABLE
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    data_[i] /= other.data_[i];
                }
            });
#elif HAS_EXECUTION
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      data_.begin(), std::divides<T>());
#else
        for (size_t i = 0; i < data_.size(); ++i) {
            data_[i] /= other.data_[i];
        }
#endif
        return *this;
    }

    /// Parallel elementwise addition of two multivectors
    /// @param other The other multivector to add
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
#if !defined(TBB_AVAILABLE) && !defined(HAS_EXECUTION) && !defined(OMP_AVAILABLE)
    [[deprecated("Consider enabling TBB, C++17 execution policies, or OpenMP for better performance.")]]
#endif
    MultiVector<T, N> parallel_add(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
#if defined(TBB_AVAILABLE)
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = data_[i] + other.data_[i];
                }
            });
#elif defined(HAS_EXECUTION)
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::plus<T>());
#elif defined(OMP_AVAILABLE)
        #pragma omp parallel for
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] + other.data_[i];
        }
#else
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::plus<T>());
#endif
        return result;
    }

    /// Parallel elementwise subtraction of two multivectors
    /// @param other The other multivector to subtract
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
#if !defined(TBB_AVAILABLE) && !defined(HAS_EXECUTION) && !defined(OMP_AVAILABLE)
    [[deprecated("Consider enabling TBB, C++17 execution policies, or OpenMP for better performance.")]]
#endif
    MultiVector<T, N> parallel_subtract(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
#if defined(TBB_AVAILABLE)
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = data_[i] - other.data_[i];
                }
            });
#elif defined(HAS_EXECUTION)
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::minus<T>());
#elif defined(OMP_AVAILABLE)
        #pragma omp parallel for
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] - other.data_[i];
        }
#else
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::minus<T>());
#endif
        return result;
    }

    /// Parallel elementwise multiplication of two multivectors
    /// @param other The other multivector to multiply with
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
#if !defined(TBB_AVAILABLE) && !defined(HAS_EXECUTION) && !defined(OMP_AVAILABLE)
    [[deprecated("Consider enabling TBB, C++17 execution policies, or OpenMP for better performance.")]]
#endif
    MultiVector<T, N> parallel_multiply(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
#if defined(TBB_AVAILABLE)
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = data_[i] * other.data_[i];
                }
            });
#elif defined(HAS_EXECUTION)
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::multiplies<T>());
#elif defined(OMP_AVAILABLE)
        #pragma omp parallel for
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] * other.data_[i];
        }
#else
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::multiplies<T>());
#endif
        return result;
    }

    /// Parallel elementwise division of two multivectors
    /// @param other The other multivector to divide by
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
#if !defined(TBB_AVAILABLE) && !defined(HAS_EXECUTION) && !defined(OMP_AVAILABLE)
    [[deprecated("Consider enabling TBB, C++17 execution policies, or OpenMP for better performance.")]]
#endif
    MultiVector<T, N> parallel_divide(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
#if defined(TBB_AVAILABLE)
        tbb::parallel_for(tbb::blocked_range<size_t>(0, data_.size()),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = data_[i] / other.data_[i];
                }
            });
#elif defined(HAS_EXECUTION)
        std::transform(std::execution::par, data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::divides<T>());
#elif defined(OMP_AVAILABLE)
        #pragma omp parallel for
        for (size_t i = 0; i < data_.size(); ++i) {
            result.data_[i] = data_[i] / other.data_[i];
        }
#else
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::divides<T>());
#endif
        return result;
    }

    /// Apply softmax transformation to each inner slice of the MultiVector
    /// @param is_log_values Whether the input values are in log space
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A new MultiVector containing the softmax transformed values
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> softmax(bool is_log_values = true, ExecutionPolicy policy = std::execution::seq) const {
        MultiVector<T, N> result(dimensions_);
        
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    if (is_log_values) {
                        // For log values, we need to compute logsumexp first
                        T max_val = *std::max_element(policy, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        T log_sum = std::log(sum) + max_val;
                        
                        // Subtract log_sum from each element
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = *it - log_sum;
                        }
                    } else {
                        // For regular values, compute softmax directly
                        T max_val = *std::max_element(policy, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        
                        // Normalize by sum
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = std::exp(*it - max_val) / sum;
                        }
                    }
                }
            }
        }
        // Optimized implementation for 2D case
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
            for (size_t i = 0; i < dim0; ++i) {
                // Calculate the start index for this row
                const size_t start_idx = i * stride0;
                
                // Get iterators for the innermost dimension
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                auto result_begin = result.data_.begin() + start_idx;
                
                if (is_log_values) {
                    // For log values, we need to compute logsumexp first
                    T max_val = *std::max_element(policy, begin, end);
                    T sum = T{};
                    for (auto it = begin; it != end; ++it) {
                        sum += std::exp(*it - max_val);
                    }
                    T log_sum = std::log(sum) + max_val;
                    
                    // Subtract log_sum from each element
                    for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                        *rit = *it - log_sum;
                    }
                } else {
                    // For regular values, compute softmax directly
                    T max_val = *std::max_element(policy, begin, end);
                    T sum = T{};
                    for (auto it = begin; it != end; ++it) {
                        sum += std::exp(*it - max_val);
                    }
                    
                    // Normalize by sum
                    for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                        *rit = std::exp(*it - max_val) / sum;
                    }
                }
            }
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto softmax_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Apply softmax for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, result_end] = result.inner_iterators(indices);
                    
                    if (is_log_values) {
                        // For log values, we need to compute logsumexp first
                        T max_val = *std::max_element(policy, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        T log_sum = std::log(sum) + max_val;
                        
                        // Subtract log_sum from each element
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = *it - log_sum;
                        }
                    } else {
                        // For regular values, compute softmax directly
                        T max_val = *std::max_element(policy, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        
                        // Normalize by sum
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = std::exp(*it - max_val) / sum;
                        }
                    }
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            softmax_recursive(softmax_recursive, indices, 0);
        }
        
        return result;
    }

    /// Apply softmax transformation to each inner slice of the MultiVector in parallel
    /// @param is_log_values Whether the input values are in log space
    /// @return A new MultiVector containing the softmax transformed values
    MultiVector<T, N> parallel_softmax(bool is_log_values = true) const {
        MultiVector<T, N> result(dimensions_);
        
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for with a 2D blocked range for better load balancing
            tbb::parallel_for(
                tbb::blocked_range2d<size_t>(0, dim0, 0, dim1),
                [&](const tbb::blocked_range2d<size_t>& range) {
                    for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
                        for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
                            // Calculate the start index for this 2D slice
                            const size_t start_idx = i * stride0 + j * stride1;
                            
                            // Get iterators for the innermost dimension
                            auto begin = data_.begin() + start_idx;
                            auto end = begin + dim2;
                            auto result_begin = result.data_.begin() + start_idx;
                            
                            if (is_log_values) {
                                // For log values, we need to compute logsumexp first
                                T max_val = *std::max_element(std::execution::seq, begin, end);
                                T sum = T{};
                                for (auto it = begin; it != end; ++it) {
                                    sum += std::exp(*it - max_val);
                                }
                                T log_sum = std::log(sum) + max_val;
                                
                                // Subtract log_sum from each element
                                for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                                    *rit = *it - log_sum;
                                }
                            } else {
                                // For regular values, compute softmax directly
                                T max_val = *std::max_element(std::execution::seq, begin, end);
                                T sum = T{};
                                for (auto it = begin; it != end; ++it) {
                                    sum += std::exp(*it - max_val);
                                }
                                
                                // Normalize by sum
                                for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                                    *rit = std::exp(*it - max_val) / sum;
                                }
                            }
                        }
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<std::pair<size_t, size_t>> indices;
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    indices.emplace_back(i, j);
                }
            }
            
            std::for_each(std::execution::par, indices.begin(), indices.end(),
                [&](const auto& idx_pair) {
                    const size_t i = idx_pair.first;
                    const size_t j = idx_pair.second;
                    
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    if (is_log_values) {
                        // For log values, we need to compute logsumexp first
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        T log_sum = std::log(sum) + max_val;
                        
                        // Subtract log_sum from each element
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = *it - log_sum;
                        }
                    } else {
                        // For regular values, compute softmax directly
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        
                        // Normalize by sum
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = std::exp(*it - max_val) / sum;
                        }
                    }
                }
            );
#endif
        }
        // Optimized implementation for 2D case
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for for 1D case
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, dim0),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i != range.end(); ++i) {
                        // Calculate the start index for this row
                        const size_t start_idx = i * stride0;
                        
                        // Get iterators for the innermost dimension
                        auto begin = data_.begin() + start_idx;
                        auto end = begin + dim1;
                        auto result_begin = result.data_.begin() + start_idx;
                        
                        if (is_log_values) {
                            // For log values, we need to compute logsumexp first
                            T max_val = *std::max_element(std::execution::seq, begin, end);
                            T sum = T{};
                            for (auto it = begin; it != end; ++it) {
                                sum += std::exp(*it - max_val);
                            }
                            T log_sum = std::log(sum) + max_val;
                            
                            // Subtract log_sum from each element
                            for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                                *rit = *it - log_sum;
                            }
                        } else {
                            // For regular values, compute softmax directly
                            T max_val = *std::max_element(std::execution::seq, begin, end);
                            T sum = T{};
                            for (auto it = begin; it != end; ++it) {
                                sum += std::exp(*it - max_val);
                            }
                            
                            // Normalize by sum
                            for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                                *rit = std::exp(*it - max_val) / sum;
                            }
                        }
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::par, indices.begin(), indices.end(),
                [&](size_t i) {
                    // Calculate the start index for this row
                    const size_t start_idx = i * stride0;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim1;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    if (is_log_values) {
                        // For log values, we need to compute logsumexp first
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        T log_sum = std::log(sum) + max_val;
                        
                        // Subtract log_sum from each element
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = *it - log_sum;
                        }
                    } else {
                        // For regular values, compute softmax directly
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        
                        // Normalize by sum
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = std::exp(*it - max_val) / sum;
                        }
                    }
                }
            );
#endif
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto softmax_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Apply softmax for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, result_end] = result.inner_iterators(indices);
                    
                    if (is_log_values) {
                        // For log values, we need to compute logsumexp first
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        T log_sum = std::log(sum) + max_val;
                        
                        // Subtract log_sum from each element
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = *it - log_sum;
                        }
                    } else {
                        // For regular values, compute softmax directly
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        
                        // Normalize by sum
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = std::exp(*it - max_val) / sum;
                        }
                    }
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            softmax_recursive(softmax_recursive, indices, 0);
        }
        
        return result;
    }

    /// Apply a unary transformation to each element of the MultiVector
    /// @param op The unary operation to apply to each element
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A new MultiVector with the transformed values
    template<typename UnaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> transform(UnaryOp op, ExecutionPolicy policy = std::execution::seq) const {
        MultiVector<T, N> result(dimensions_);
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    // Apply transformation to each element
                    std::transform(policy, begin, end, result_begin, op);
                }
            }
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
            for (size_t i = 0; i < dim0; ++i) {
                // Calculate the start index for this row
                const size_t start_idx = i * stride0;
                
                // Get iterators for the innermost dimension
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                auto result_begin = result.data_.begin() + start_idx;
                
                // Apply transformation to each element
                std::transform(policy, begin, end, result_begin, op);
            }
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Apply transformation for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    std::transform(policy, begin, end, result_begin, op);
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            transform_recursive(transform_recursive, indices, 0);
        }
        return result;
    }

    /// Apply a unary transformation to each element of the MultiVector in parallel
    /// @param op The unary operation to apply to each element
    /// @return A new MultiVector with the transformed values
    template<typename UnaryOp>
    MultiVector<T, N> parallel_transform(UnaryOp op) const {
        MultiVector<T, N> result(dimensions_);
        // Optimized implementation for 3D case (most common)
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            
            // Pre-calculate strides for faster access
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for with a 2D blocked range for better load balancing
            tbb::parallel_for(
                tbb::blocked_range2d<size_t>(0, dim0, 0, dim1),
                [&](const tbb::blocked_range2d<size_t>& range) {
                    for (size_t i = range.rows().begin(); i != range.rows().end(); ++i) {
                        for (size_t j = range.cols().begin(); j != range.cols().end(); ++j) {
                            // Calculate the start index for this 2D slice
                            const size_t start_idx = i * stride0 + j * stride1;
                            
                            // Get iterators for the innermost dimension
                            auto begin = data_.begin() + start_idx;
                            auto end = begin + dim2;
                            auto result_begin = result.data_.begin() + start_idx;
                            
                            // Apply transformation to each element
                            std::transform(std::execution::seq, begin, end, result_begin, op);
                        }
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0 * dim1);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::seq, 
                indices.begin(), 
                indices.end(),
                [&](size_t idx) {
                    const size_t i = idx / dim1;
                    const size_t j = idx % dim1;
                    
                    // Calculate the start index for this 2D slice
                    const size_t start_idx = i * stride0 + j * stride1;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    // Apply transformation to each element
                    std::transform(std::execution::seq, begin, end, result_begin, op);
                }
            );
#endif
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
#if TBB_AVAILABLE
            // Use TBB's parallel_for for 1D case
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, dim0),
                [&](const tbb::blocked_range<size_t>& range) {
                    for (size_t i = range.begin(); i != range.end(); ++i) {
                        // Calculate the start index for this row
                        const size_t start_idx = i * stride0;
                        
                        // Get iterators for the innermost dimension
                        auto begin = data_.begin() + start_idx;
                        auto end = begin + dim1;
                        auto result_begin = result.data_.begin() + start_idx;
                        
                        // Apply transformation to each element
                        std::transform(std::execution::seq, begin, end, result_begin, op);
                    }
                }
            );
#else
            // Fallback to standard library parallel execution
            std::vector<size_t> indices(dim0);
            std::iota(indices.begin(), indices.end(), 0);
            std::for_each(std::execution::seq, 
                indices.begin(), 
                indices.end(),
                [&](size_t i) {
                    // Calculate the start index for this row
                    const size_t start_idx = i * stride0;
                    
                    // Get iterators for the innermost dimension
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim1;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    // Apply transformation to each element
                    std::transform(std::execution::seq, begin, end, result_begin, op);
                }
            );
#endif
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Apply transformation for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    std::transform(std::execution::seq, begin, end, result_begin, op);
                    return;
                }
                
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    indices[dim] = i;
                    self(self, indices, dim + 1);
                }
            };

            std::array<size_t, N-1> indices{};
            transform_recursive(transform_recursive, indices, 0);
        }
        return result;
    }

    // Convenience methods for common unary operations
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> negate(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return -x; }, policy);
    }

    MultiVector<T, N> parallel_negate() const {
        return parallel_transform([](const T& x) { return -x; });
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> abs(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return std::abs(x); }, policy);
    }

    MultiVector<T, N> parallel_abs() const {
        return parallel_transform([](const T& x) { return std::abs(x); });
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> exp(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return std::exp(x); }, policy);
    }

    MultiVector<T, N> parallel_exp() const {
        return parallel_transform([](const T& x) { return std::exp(x); });
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> log(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return std::log(x); }, policy);
    }

    MultiVector<T, N> parallel_log() const {
        return parallel_transform([](const T& x) { return std::log(x); });
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, N> sqrt(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return std::sqrt(x); }, policy);
    }

    MultiVector<T, N> parallel_sqrt() const {
        return parallel_transform([](const T& x) { return std::sqrt(x); });
    }

private:
    std::vector<T> data_;
    std::array<size_t, N> dimensions_;
    std::array<size_t, N> strides_;

    // Helper function to calculate the 1D index from multi-dimensional indices
    /// @param indices The multi-dimensional indices.
    /// @return The calculated 1D index.
    /// @throws std::out_of_range if any index is out of bounds.
    inline size_t calculate_index(const std::array<size_t, N>& indices) const {
#ifndef NDEBUG
        for (size_t i = 0; i < indices.size(); ++i) {
            if (indices[i] >= dimensions_[i]) {
                throw std::out_of_range("Index out of bounds");
            }
        }
#endif
        // Optimized implementation for 3D case
        if constexpr (N == 3) {
            return indices[0] * strides_[0] + indices[1] * strides_[1] + indices[2];
        } else if constexpr (N == 2) {
            return indices[0] * strides_[0] + indices[1];
        } else if constexpr (N == 1) {
            return indices[0];
        } else {
            size_t index = 0;
            const size_t* stride = strides_.data();
            const size_t* idx = indices.data();
            for (size_t i = N; i-- > 0;) {
                index += idx[i] * stride[i];
            }
            return index;
        }
    }

    /// Helper function to calculate start and end indices for iterators
    /// @param outer_indices The indices of the outer dimensions.
    /// @return A pair of start and end indices for the iterators.
    std::pair<size_t, size_t> calculate_start_end_indices(const std::array<size_t, N - 1>& outer_indices) const {
        if constexpr (N == 1) {
            return {0, dimensions_.back()};
        } else {
            // Optimized implementation for 3D case
            if constexpr (N == 3) {
                const size_t start_index = outer_indices[0] * strides_[0] + outer_indices[1] * strides_[1];
                const size_t end_index = start_index + dimensions_.back();
                return {start_index, end_index};
            } else {
                std::array<size_t, N> full_indices{0};
                std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
                const size_t start_index = calculate_index(full_indices);
                const size_t end_index = start_index + dimensions_.back();
                return {start_index, end_index};
            }
        }
    }
};

/// Specialization of MultiVector for N=1 that returns a scalar from reduce functions
template <typename T>
class MultiVector<T, 1> {
public:
    /// Constructor that takes dimensions
    /// @param dimensions The dimensions of the MultiVector.
    /// @throws std::invalid_argument if dimensions are empty.
    MultiVector(const std::array<size_t, 1>& dimensions) : dimensions_(dimensions) {
#ifndef NDEBUG
        if (dimensions.empty()) {
            throw std::invalid_argument("Dimensions cannot be empty");
        }
#endif
        // Calculate total size and strides
        strides_[0] = 1;
        data_.resize(dimensions_[0], T{});
    }

    MultiVector() {
        dimensions_ = {0};
        data_ = {};
        strides_ = {1};
    }

    void resize(const std::array<size_t, 1>& dimensions) {
#ifndef NDEBUG
        if (dimensions.empty()) {
            throw std::invalid_argument("Dimensions cannot be empty");
        }
#endif
        dimensions_ = dimensions;
        strides_[0] = 1;
        data_.resize(dimensions_[0], T{});
    }

    void resize(const std::array<size_t, 1>& dimensions, const T& value) {
#ifndef NDEBUG
        if (dimensions.empty()) {
            throw std::invalid_argument("Dimensions cannot be empty");
        }
#endif
        dimensions_ = dimensions;
        strides_[0] = 1;
        data_.resize(dimensions_[0], value);
    }

    /// Get a span view of the data
    /// @return A span view of the data
    std::span<T> as_span() {
        return std::span<T>(data_);
    }

    /// Get a const span view of the data
    /// @return A const span view of the data
    std::span<const T> as_span() const {
        return std::span<const T>(data_);
    }

    /// Implicit conversion to a span
    /// @return A span view of the data
    operator std::span<T>() {
        return as_span();
    }

    /// Implicit conversion to a const span
    /// @return A const span view of the data
    operator std::span<const T>() const {
        return as_span();
    }
    
    /// Fill the vector with a single value
    /// @param value The value to fill the vector with.
    void inner_fill(const T& value) {
        std::fill(data_.begin(), data_.end(), value);
    }

    /// Fill the vector with a range of values
    /// @param values The values to fill the vector with.
    void inner_fill(const std::span<T const> values) {
#ifndef NDEBUG
        if (values.size() > dimensions_[0]) {
            throw std::invalid_argument("(" + std::to_string(1) + "D) Incorrect number of values, expected " + std::to_string(dimensions_[0]) + " but got " + std::to_string(values.size()));
        }
#endif
        std::copy(values.begin(), values.end(), data_.begin());
    }

    /// Access element at a given index
    /// @param indices The index of the element to access.
    /// @return A reference to the element at the specified index.
    /// @throws std::invalid_argument if the number of indices is incorrect.
    T& at(const std::array<size_t, 1>& indices) {
#ifndef NDEBUG
        if (indices.size() != dimensions_.size()) {
            throw std::invalid_argument("Incorrect number of indices");
        }
#endif
        return data_.at(calculate_index(indices));
    }

    const T& at(const std::array<size_t, 1>& indices) const {
#ifndef NDEBUG
        if (indices.size() != dimensions_.size()) {
            throw std::invalid_argument("Incorrect number of indices");
        }
#endif
        return data_.at(calculate_index(indices));
    }

    /// Get the size of the dimension
    /// @return The size of the dimension.
    size_t size(size_t dimension) const {
#ifndef NDEBUG
        if (dimension >= dimensions_.size()) {
            throw std::out_of_range("Dimension out of range");
        }
#endif
        return dimensions_[dimension];
    }

    /// Get the total number of elements
    /// @return The total number of elements in the MultiVector.
    size_t total_size() const {
        return data_.size();
    }

    /// Iterator for the vector
    /// @return An iterator to the beginning of the vector.
    typename std::vector<T>::iterator inner_begin() {
        return data_.begin();
    }

    typename std::vector<T>::iterator inner_end() {
        return data_.end();
    }

    /// Const iterator for the vector
    /// @return A const iterator to the beginning of the vector.
    typename std::vector<T>::const_iterator inner_begin() const {
        return data_.cbegin();
    }

    typename std::vector<T>::const_iterator inner_end() const {
        return data_.cend();
    }

    std::pair<typename std::vector<T>::const_iterator, typename std::vector<T>::const_iterator> inner_iterators() const {
        return {data_.cbegin(), data_.cend()};
    }

    /// Clear the MultiVector
    /// @note This will clear the data contained in the MultiVector.
    void clear() {
        data_.clear();
    }

    std::array<size_t, 1> dimensions() const {
        return dimensions_;
    }

    const std::array<size_t, 1>& strides() const {
        return strides_;
    }

    const std::vector<T>& data() const {
        return data_;
    }

    /// Reduce the vector using a binary operation, returning a scalar
    /// @param binary_op The binary operation to use for reduction
    /// @param init The initial value for the reduction
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A scalar value containing the reduced value
    template<typename BinaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    auto reduce(BinaryOp binary_op, const auto& init, ExecutionPolicy policy = std::execution::seq) const {
        if constexpr (std::is_same_v<ExecutionPolicy, std::execution::sequenced_policy>) {
            return std::accumulate(data_.begin(), data_.end(), init, binary_op);
        } else {
            return std::reduce(policy, data_.begin(), data_.end(), init, binary_op);
        }
    }

    /// Optimized parallel reduce using TBB for better performance if available, otherwise falls back to standard library
    /// @param binary_op The binary operation to use for reduction
    /// @param init The initial value for the reduction
    /// @return A scalar value containing the reduced value
    template<typename BinaryOp>
    auto parallel_reduce(BinaryOp binary_op, const auto& init) const {
        using ResultType = decltype(binary_op(init, std::declval<T>()));
#if TBB_AVAILABLE
        // Use TBB's parallel_reduce for better performance
        return tbb::parallel_reduce(
            tbb::blocked_range<typename std::vector<T>::const_iterator>(data_.begin(), data_.end()),
            init,
            [&binary_op](const tbb::blocked_range<typename std::vector<T>::const_iterator>& range, ResultType init) {
                return std::reduce(std::execution::seq, range.begin(), range.end(), init, binary_op);
            },
            binary_op
        );
#else
        // Fallback to standard library execution
        return std::reduce(std::execution::seq, data_.begin(), data_.end(), init, binary_op);
#endif
    }

    // Convenience methods for common reductions with optional execution policy
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T sum(ExecutionPolicy policy = std::execution::seq) const {
        return reduce(std::plus<T>(), T{}, policy);
    }

    // Optimized parallel sum method
    T parallel_sum() const {
        return parallel_reduce(std::plus<T>(), T{});
    }

    T full_sum() const {
        return sum();
    }

    T parallel_full_sum() const {
        return parallel_sum();
    }

    T full_product() const {
        return product();
    }

    T parallel_full_product() const {
        return parallel_product();
    }

    T full_min() const {
        return min();
    }

    T parallel_full_min() const {
        return parallel_min();
    }

    T full_max() const {
        return max();
    }

    T parallel_full_max() const {
        return parallel_max();
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T product(ExecutionPolicy policy = std::execution::seq) const {
        return reduce(std::multiplies<T>(), T{1}, policy);
    }

    // Optimized parallel product method
    T parallel_product() const {
        return parallel_reduce(std::multiplies<T>(), T{1});
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T max(ExecutionPolicy policy = std::execution::seq) const {
        return reduce([](const T& a, const T& b) { return std::max(a, b); },
                     std::numeric_limits<T>::lowest(), policy);
    }

    // Optimized parallel max method
    T parallel_max() const {
        return parallel_reduce([](const T& a, const T& b) { return std::max(a, b); },
                              std::numeric_limits<T>::lowest());
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T min(ExecutionPolicy policy = std::execution::seq) const {
        return reduce([](const T& a, const T& b) { return std::min(a, b); },
                     std::numeric_limits<T>::max(), policy);
    }

    // Optimized parallel min method
    T parallel_min() const {
        return parallel_reduce([](const T& a, const T& b) { return std::min(a, b); },
                              std::numeric_limits<T>::max());
    }

    /// Compute the log-sum-exp reduction over the vector
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A scalar value containing the log-sum-exp value
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T logsumexp(ExecutionPolicy policy = std::execution::seq) const {
        T max_val = *std::max_element(policy, data_.begin(), data_.end());
        T sum = std::transform_reduce(
            policy,
            data_.begin(), data_.end(),
            T{},
            std::plus<T>(),
            [max_val](const T& x) { return std::exp(x - max_val); }
        );
        return std::log(sum) + max_val;
    }

    /// Compute the log-sum-exp reduction over the vector in parallel
    /// @return A scalar value containing the log-sum-exp value
    T parallel_logsumexp() const {
        #if TBB_AVAILABLE
        auto policy = std::execution::par;
        #else
        auto policy = std::execution::seq;
        #endif

        T max_val = *std::max_element(policy, data_.begin(), data_.end());
        T sum = std::transform_reduce(
            policy,
            data_.begin(), data_.end(),
            T{},
            std::plus<T>(),
            [max_val](const T& x) { return std::exp(x - max_val); }
        );
        return std::log(sum) + max_val;
    }

    /// Apply a transformation to each element of the vector
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element
    /// @param execution_policy Optional execution policy for parallel execution
    template<typename BinaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    void element_transform(BinaryOp op, const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) {
#ifndef NDEBUG
        if (values.size() != dimensions_[0]) {
            throw std::invalid_argument("Values span size does not match the vector size");
        }
#endif
        std::transform(policy, data_.begin(), data_.end(), values.begin(), data_.begin(), op);
    }

    /// Apply a transformation to each element of the vector in parallel
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element
    template<typename BinaryOp>
    void parallel_element_transform(BinaryOp op, const std::span<T const> values) {
#ifndef NDEBUG
        if (values.size() != dimensions_[0]) {
            throw std::invalid_argument("Values span size does not match the vector size");
        }
#endif
#if TBB_AVAILABLE
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, dimensions_[0]),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    data_[i] = op(data_[i], values[i]);
                }
            }
        );
#else
        std::transform(std::execution::seq, data_.begin(), data_.end(), values.begin(), data_.begin(), op);
#endif
    }

    // Convenience methods for element-level operations
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> element_add(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::plus<T>(), values, policy);
    }

    MultiVector<T, 1> parallel_element_add(const std::span<T const> values) const {
        return parallel_transform(std::plus<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> element_subtract(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::minus<T>(), values, policy);
    }

    MultiVector<T, 1> parallel_element_subtract(const std::span<T const> values) const {
        return parallel_transform(std::minus<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> element_multiply(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::multiplies<T>(), values, policy);
    }

    MultiVector<T, 1> parallel_element_multiply(const std::span<T const> values) const {
        return parallel_transform(std::multiplies<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> element_divide(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::divides<T>(), values, policy);
    }

    MultiVector<T, 1> parallel_element_divide(const std::span<T const> values) const {
        return parallel_transform(std::divides<T>(), values);
    }

    std::vector<T> data_;
    std::array<size_t, 1> dimensions_;
    std::array<size_t, 1> strides_;

    // Helper function to calculate the 1D index from multi-dimensional indices
    /// @param indices The multi-dimensional indices.
    /// @return The calculated 1D index.
    /// @throws std::out_of_range if any index is out of bounds.
    inline size_t calculate_index(const std::array<size_t, 1>& indices) const {
#ifndef NDEBUG
        if (indices[0] >= dimensions_[0]) {
            throw std::out_of_range("Index out of bounds");
        }
#endif
        return indices[0];
    }

    /// Helper function to calculate start and end indices for iterators
    /// @param outer_indices The indices of the outer dimensions.
    /// @return A pair of start and end indices for the iterators.
    std::pair<size_t, size_t> calculate_start_end_indices(const std::array<size_t, 0>& outer_indices) const {
        return {0, dimensions_[0]};
    }

    /// Apply a unary transformation to each element of the vector
    /// @param op The unary operation to apply to each element
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A new MultiVector with the transformed values
    template<typename UnaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> transform(UnaryOp op, ExecutionPolicy policy = std::execution::seq) const {
        MultiVector<T, 1> result(dimensions_);
        std::transform(policy, data_.begin(), data_.end(), result.data_.begin(), op);
        return result;
    }

    /// Apply a unary transformation to each element of the vector in parallel
    /// @param op The unary operation to apply to each element
    /// @return A new MultiVector with the transformed values
    template<typename UnaryOp>
    MultiVector<T, 1> parallel_transform(UnaryOp op) const {
        MultiVector<T, 1> result(dimensions_);
#if TBB_AVAILABLE
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, dimensions_[0]),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = op(data_[i]);
                }
            }
        );
#else
        std::transform(std::execution::seq, data_.begin(), data_.end(), result.data_.begin(), op);
#endif
        return result;
    }

    // Convenience methods for common unary operations
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> negate(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return -x; }, policy);
    }

    MultiVector<T, 1> parallel_negate() const {
        return parallel_transform([](const T& x) { return -x; });
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> abs(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return std::abs(x); }, policy);
    }

    MultiVector<T, 1> parallel_abs() const {
        return parallel_transform([](const T& x) { return std::abs(x); });
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> exp(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return std::exp(x); }, policy);
    }

    MultiVector<T, 1> parallel_exp() const {
        return parallel_transform([](const T& x) { return std::exp(x); });
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> log(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return std::log(x); }, policy);
    }

    MultiVector<T, 1> parallel_log() const {
        return parallel_transform([](const T& x) { return std::log(x); });
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> sqrt(ExecutionPolicy policy = std::execution::seq) const {
        return transform([](const T& x) { return std::sqrt(x); }, policy);
    }

    MultiVector<T, 1> parallel_sqrt() const {
        return parallel_transform([](const T& x) { return std::sqrt(x); });
    }

    /// Apply a transformation to each element of the vector
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A new MultiVector with the transformed values
    template<typename BinaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> transform(BinaryOp op, const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
#ifndef NDEBUG
        if (values.size() != dimensions_[0]) {
            throw std::invalid_argument("Values span size does not match the vector size");
        }
#endif
        MultiVector<T, 1> result(dimensions_);
        std::transform(policy, data_.begin(), data_.end(), values.begin(), result.data_.begin(), op);
        return result;
    }

    /// Apply a transformation to each element of the vector in parallel
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element
    /// @return A new MultiVector with the transformed values
    template<typename BinaryOp>
    MultiVector<T, 1> parallel_transform(BinaryOp op, const std::span<T const> values) const {
#ifndef NDEBUG
        if (values.size() != dimensions_[0]) {
            throw std::invalid_argument("Values span size does not match the vector size");
        }
#endif
        MultiVector<T, 1> result(dimensions_);
#if TBB_AVAILABLE
        tbb::parallel_for(
            tbb::blocked_range<size_t>(0, dimensions_[0]),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    result.data_[i] = op(data_[i], values[i]);
                }
            }
        );
#else
        std::transform(std::execution::seq, data_.begin(), data_.end(), values.begin(), result.data_.begin(), op);
#endif
        return result;
    }

    // Convenience methods for element-level operations
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> add(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::plus<T>(), values, policy);
    }

    MultiVector<T, 1> parallel_add(const std::span<T const> values) const {
        return parallel_transform(std::plus<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> subtract(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::minus<T>(), values, policy);
    }

    MultiVector<T, 1> parallel_subtract(const std::span<T const> values) const {
        return parallel_transform(std::minus<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> multiply(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::multiplies<T>(), values, policy);
    }

    MultiVector<T, 1> parallel_multiply(const std::span<T const> values) const {
        return parallel_transform(std::multiplies<T>(), values);
    }

    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    MultiVector<T, 1> divide(const std::span<T const> values, ExecutionPolicy policy = std::execution::seq) const {
        return transform(std::divides<T>(), values, policy);
    }

    MultiVector<T, 1> parallel_divide(const std::span<T const> values) const {
        return parallel_transform(std::divides<T>(), values);
    }

};


template <typename T, size_t N>
/// RaggedMultiVector class template
/// This class represents a multi-dimensional vector with a ragged internal dimension.
/// It provides methods to access elements, get sizes, and iterate over dimensions.
class RaggedMultiVector {
public:
    /// Constructor that takes dimensions and ragged dimensions
    /// @param dimensions The dimensions of the MultiVector.
    /// @param ragged_dimensions The ragged dimensions of the MultiVector.
    /// @param default_value The default value for elements.
    RaggedMultiVector(const std::array<size_t, N - 1>& dimensions, const std::span<size_t const> ragged_dimensions, const T& default_value = T{}) : dimensions_(dimensions), ragged_dimensions_(ragged_dimensions.begin(), ragged_dimensions.end()) {
        // Calculate total size and strides
        if constexpr (N == 1) {
            throw std::invalid_argument("RaggedMultiVector does not support 1D vectors");
        } else {
            size_t total_size = 0;
            size_t data_size = 0;
            size_t total_entries = 1;

            std::vector<size_t> entry_dimensions(dimensions_.begin(), dimensions_.end() - 1);

            // Calculate the product of the first N - 1 dimensions corresponding to the total number of entries
            total_entries = std::reduce(std::execution::seq, entry_dimensions.begin(), entry_dimensions.end(), 1, std::multiplies<size_t>());

            data_size = std::reduce(std::execution::seq, ragged_dimensions_.begin(), ragged_dimensions_.end());

            total_size = total_entries * data_size;

            // Calculate the strides for the N - 2 dimensions
            strides_[N - 3] = data_size;
            for (size_t i = N - 3; i > 0; i--) {
                strides_[i - 1] = strides_[i] * entry_dimensions[i];
            }

            // ragged strides are the cumulative sum of the ragged dimensions and provide the offset necessary to access the ragged dimension.
            ragged_offsets_.resize(ragged_dimensions_.size(), 0);
            for (size_t i = 1; i < ragged_dimensions_.size(); ++i) {
                ragged_offsets_[i] = ragged_offsets_[i - 1] + ragged_dimensions_[i - 1];
            }

            data_.resize(total_size, default_value);
        }
    }


    RaggedMultiVector() {
        dimensions_ = {0};
        ragged_dimensions_ = {0};
        ragged_offsets_ = {0};
        data_ = {0};
    }

    void resize(const std::array<size_t, N - 1>& dimensions, const std::span<size_t const> ragged_dimensions, const T& default_value = T{}) {
        dimensions_ = dimensions;
        ragged_dimensions_ = std::vector<size_t>(ragged_dimensions.begin(), ragged_dimensions.end());

        size_t total_size = 0;
        size_t data_size = 0;
        size_t total_entries = 1;

        std::vector<size_t> entry_dimensions(dimensions_.begin(), dimensions_.end() - 1);

        total_entries = std::reduce(std::execution::seq, entry_dimensions.begin(), entry_dimensions.end(), 1, std::multiplies<size_t>());

        data_size = std::reduce(std::execution::seq, ragged_dimensions_.begin(), ragged_dimensions_.end());

        total_size = total_entries * data_size;

        // Calculate the strides for the N - 2 dimensions
        strides_[N - 3] = data_size;
        for (size_t i = N - 3; i > 0; i--) {
            strides_[i - 1] = strides_[i] * entry_dimensions[i];
        }

        ragged_offsets_.resize(ragged_dimensions_.size(), 0);
        for (size_t i = 1; i < ragged_dimensions_.size(); ++i) {
            ragged_offsets_[i] = ragged_offsets_[i - 1] + ragged_dimensions_[i - 1];
        }

        data_.resize(total_size, default_value);
    }


    void inner_fill(const std::array<size_t, N - 1>& indices, const T& value) {
        const auto [begin, end] = inner_iterators(indices);
        std::fill(begin, end, value);
    }

    void inner_fill(const std::array<size_t, N - 1>& indices, const std::span<T const> values) {
#ifndef NDEBUG
        if (values.size() > ragged_dimensions_.at(indices.back())) {
            throw std::invalid_argument("(" + std::to_string(N) + "D) Incorrect number of values, expected " + std::to_string(ragged_dimensions_.at(indices.back())) + " but got " + std::to_string(values.size()));
        }
#endif
        std::copy(values.begin(), values.end(), inner_begin(indices));
    }

    size_t total_size() const {
        return data_.size();
    }

    /// Access element at a given multi-dimensional index
    /// @param indices The indices of the element to access.
    /// @return A reference to the element at the specified indices.
    T& at(const std::array<size_t, N>& indices) {
#ifndef NDEBUG
        if (indices.size() != N) {
            throw std::invalid_argument("Incorrect number of indices");
        }
#endif
        return data_.at(calculate_index(indices));
    }

    const T& at(const std::array<size_t, N>& indices) const {
#ifndef NDEBUG
        if (indices.size() != N) {
            throw std::invalid_argument("Incorrect number of indices");
        }
#endif
        return data_.at(calculate_index(indices));
    }

    /// Iterator for the innermost dimension
    /// @param outer_indices The indices of the outer dimensions.
    /// @return An iterator to the beginning of the innermost dimension.
    typename std::vector<T>::iterator inner_begin(const std::array<size_t, N - 1>& outer_indices) {
        auto [start_index, _] = calculate_start_end_indices(outer_indices);
        return data_.begin() + start_index;
    }

    typename std::vector<T>::iterator inner_end(const std::array<size_t, N - 1>& outer_indices) {
        auto [_, end_index] = calculate_start_end_indices(outer_indices);
        return data_.begin() + end_index;
    }

    std::pair<typename std::vector<T>::iterator, typename std::vector<T>::iterator> inner_iterators(const std::array<size_t, N - 1>& outer_indices) {
        auto [start_index, end_index] = calculate_start_end_indices(outer_indices);
        return {data_.begin() + start_index, data_.begin() + end_index};
    }

    /// Const iterator for the innermost dimension
    /// @param outer_indices The indices of the outer dimensions.
    /// @return A const iterator to the beginning of the innermost dimension.
    /// @throws std::invalid_argument if the number of outer indices is incorrect.
    typename std::vector<T>::const_iterator inner_begin(const std::array<size_t, N - 1>& outer_indices) const {
        auto [start_index, _] = calculate_start_end_indices(outer_indices);
        return data_.cbegin() + start_index;
    }

    typename std::vector<T>::const_iterator inner_end(const std::array<size_t, N - 1>& outer_indices) const {
        auto [_, end_index] = calculate_start_end_indices(outer_indices);
        return data_.cbegin() + end_index;
    }

    std::pair<typename std::vector<T>::const_iterator, typename std::vector<T>::const_iterator> inner_iterators(const std::array<size_t, N - 1>& outer_indices) const {
        auto [start_index, end_index] = calculate_start_end_indices(outer_indices);
        return {data_.cbegin() + start_index, data_.cbegin() + end_index};
    }

    void clear() {
        data_.clear();
    }

    std::array<size_t, N> dimensions() const {
        std::array<size_t, N> dimensions;
        std::copy(dimensions_.begin(), dimensions_.end(), dimensions.begin());
        dimensions.back() = ragged_dimensions_.size();
        return dimensions;
    }

    std::vector<size_t> ragged_dimensions() const {
        return ragged_dimensions_;
    }

    size_t ragged_dimensions(size_t index) const {
        return ragged_dimensions_.at(index);
    }

    const std::vector<T>& data() const {
        return data_;
    }

// private:
    std::vector<T> data_;
    std::array<size_t, N - 1> dimensions_;
    std::vector<size_t> ragged_dimensions_;
    std::vector<size_t> ragged_offsets_;
    std::array<size_t, N - 2> strides_;


    /// Helper function to calculate the 1D index from multi-dimensional indices
    /// @param indices The multi-dimensional indices.
    /// @return The calculated 1D index.
    /// @throws std::out_of_range if any index is out of bounds.
    inline size_t calculate_index(const std::array<size_t, N>& indices) const {
#ifndef NDEBUG
        if constexpr (N == 1) {
            if (indices[0] >= ragged_dimensions_.size()) {
                throw std::out_of_range("Index out of bounds");
            }
        } else {
            for (size_t i = 0; i < N - 1; ++i) {
                if (indices[i] >= dimensions_[i]) {
                    throw std::out_of_range("Index out of bounds");
                }
            }
            if (indices.back() >= ragged_dimensions_.at(indices[N - 2])) {
                throw std::out_of_range("Index out of bounds");
            }
        }
#endif

        size_t index = 0;
        for (size_t i = 0; i < N - 2; ++i) {
            index += indices[i] * strides_[i];
        }
        index += ragged_offsets_[indices[N - 2]] + indices.back();
        return index;
    }

    /// Helper function to calculate start and end indices for iterators
    /// @param outer_indices The indices of the outer dimensions.
    /// @return A pair of start and end indices for the iterators.
    std::pair<size_t, size_t> calculate_start_end_indices(const std::array<size_t, N - 1>& outer_indices) const {
        if constexpr (N == 1) {
            return {0, ragged_dimensions_.size()};
        } else {
            std::array<size_t, N> full_indices{};
            std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
            const size_t start_index = calculate_index(full_indices);
            const size_t end_index = start_index + ragged_dimensions_.at(outer_indices[N - 2]);
            return {start_index, end_index};
        }
    }
};

