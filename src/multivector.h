#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <array>
#include <span>
#include <cmath>
#include <ranges>
#include <utility>

#include "multivector_ops.h"
#include "multivector_dispatch.h"
#include "multivector_algorithms.h"
#include "profiler.h"

/// MultiVector class template
/// This class represents a multi-dimensional vector with fixed dimensions.
/// It provides methods to access elements, get sizes, and iterate over dimensions.
///
/// Invariants (always maintained after construction/resize):
/// - dimensions_.size() == N; all dimensions_[i] >= 0.
/// - total_size() == product of dimensions_[i] (0 if any dimension is 0).
/// - data_.size() == total_size().
/// - strides_[i] = product of dimensions_[i+1..N-1] (strides_[N-1] == 1).
/// - For valid indices, data_[calculate_index(indices)] is the element at that position.
template <typename T, size_t N>
class MultiVector {
    static_assert(N >= 1, "MultiVector requires N >= 1");
public:
    /// Constructor that takes dimensions
    /// @param dimensions The dimensions of the MultiVector.
    MultiVector(const std::array<size_t, N>& dimensions) : dimensions_(dimensions) {
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
    /// @throws std::invalid_argument if data size does not match product of dimensions.
    MultiVector(const std::vector<T>& data, const std::array<size_t, N>& dimensions) : dimensions_(dimensions) {
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

    /// Default constructor: creates an empty (valid) MultiVector with total_size() == 0.
    /// Invariants: dimensions_ zero-initialized, strides_ and data_ consistent.
    MultiVector() {
        dimensions_ = {};
        size_t total_size = 1;
        for (size_t i = N; i-- > 0;) {
            strides_[i] = total_size;
            total_size *= dimensions_[i];
        }
        data_.resize(total_size, T{});
    }

    void resize(const std::array<size_t, N>& dimensions) {
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
        if (values.size() > dimensions_.back()) {
            throw std::invalid_argument("(" + std::to_string(N) + "D) Incorrect number of values, expected " + std::to_string(dimensions_.back()) + " but got " + std::to_string(values.size()));
        }
        auto [begin, end] = inner_iterators(indices);
        std::copy(values.begin(), values.end(), begin);
    }

    /// Access element at a given multi-dimensional index (always bounds-checked).
    /// @param indices The indices of the element to access.
    /// @return A reference to the element at the specified indices.
    /// @throws std::out_of_range if any index is out of bounds.
    T& at(const std::array<size_t, N>& indices) {
        for (size_t i = 0; i < N; ++i) {
            if (indices[i] >= dimensions_[i]) {
                throw std::out_of_range("MultiVector::at index out of bounds");
            }
        }
        return data_.at(calculate_index(indices));
    }

    const T& at(const std::array<size_t, N>& indices) const {
        for (size_t i = 0; i < N; ++i) {
            if (indices[i] >= dimensions_[i]) {
                throw std::out_of_range("MultiVector::at index out of bounds");
            }
        }
        return data_.at(calculate_index(indices));
    }

    /// Get the size of a specific dimension
    /// @param dimension The dimension to query (must be < N).
    /// @return The size of the specified dimension.
    /// @throws std::out_of_range if dimension >= N.
    size_t size(size_t dimension) const {
        if (dimension >= N) {
            throw std::out_of_range("Dimension out of range");
        }
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
    typename std::vector<T>::iterator inner_begin(const std::array<size_t, N-1>& outer_indices) {
        std::array<size_t, N> full_indices{0};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        size_t start_index = calculate_index(full_indices);
        return data_.begin() + start_index;
    }

    typename std::vector<T>::iterator inner_end(const std::array<size_t, N-1>& outer_indices) {
        std::array<size_t, N> full_indices{0};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        size_t start_index = calculate_index(full_indices);
        return data_.begin() + start_index + dimensions_.back();
    }

    /// Const iterator for the innermost dimension
    /// @param outer_indices The indices of the outer dimensions.
    /// @return A const iterator to the beginning of the innermost dimension.
    typename std::vector<T>::const_iterator inner_begin(const std::array<size_t, N-1>& outer_indices) const {
        std::array<size_t, N> full_indices{0};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        size_t start_index = calculate_index(full_indices);
        return data_.cbegin() + start_index;
    }

    typename std::vector<T>::const_iterator inner_end(const std::array<size_t, N-1>& outer_indices) const {
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

    /// Fill the MultiVector with a single value
    /// @param value The value to fill the MultiVector with.
    void fill(const T& value) {
        std::fill(data_.begin(), data_.end(), value);
    }

    /// Reduce the innermost dimension using a binary operation, returning a MultiVector of dimension N-1
    /// @param binary_op The binary operation to use for reduction
    /// @param init The initial value for the reduction
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A MultiVector<T, N-1> containing the reduced values
    template<typename BinaryOp, typename ExecutionPolicy = std::execution::sequenced_policy>
    auto reduce(BinaryOp binary_op, const auto& init, ExecutionPolicy policy = std::execution::seq) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::reduce");
#endif
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
            
            for (size_t i = 0; i < dim0; ++i) {
                for (size_t j = 0; j < dim1; ++j) {
                    const size_t start_idx = i * stride0 + j * stride1;
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    result.unchecked_at({i, j}) = moire_kernels::reduce_slice(begin, end, init, binary_op, policy);
                }
            }
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            for (size_t i = 0; i < dim0; ++i) {
                const size_t start_idx = i * stride0;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                result.unchecked_at({i}) = moire_kernels::reduce_slice(begin, end, init, binary_op, policy);
            }
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto reduce_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    // Perform reduction for these outer indices
                    const auto [begin, end] = inner_iterators(indices);
                    result.unchecked_at(indices) = moire_kernels::reduce_slice(begin, end, init, binary_op, policy);
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

    /// Optimized parallel reduce using parallel execution when enabled and workload exceeds threshold
    /// @param binary_op The binary operation to use for reduction
    /// @param init The initial value for the reduction
    /// @return A MultiVector<T, N-1> containing the reduced values
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    template<typename BinaryOp>
    auto parallel_reduce(BinaryOp binary_op, const auto& init) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::parallel_reduce");
#endif
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
            
            moire_dispatch::for_each_slice_2d(dim0, dim1, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i, size_t j) {
                const size_t start_idx = i * stride0 + j * stride1;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim2;
                result.unchecked_at({i, j}) = moire_kernels::reduce_slice(begin, end, init, binary_op, std::execution::seq);
            });
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            moire_dispatch::for_each_slice_1d(dim0, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i) {
                const size_t start_idx = i * stride0;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                result.unchecked_at({i}) = moire_kernels::reduce_slice(begin, end, init, binary_op, std::execution::seq);
            });
        }
        else {
            // Fallback to recursive approach for other dimensions
            // Helper to generate all possible outer indices
            auto reduce_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    const auto [begin, end] = inner_iterators(indices);
                    result.unchecked_at(indices) = moire_kernels::reduce_slice(begin, end, init, binary_op, std::execution::seq);
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
     * @brief Fully reduces the multivector down to a scalar value with parallel execution when enabled
     * 
     * This function applies a binary operation to all elements in the multivector,
     * reducing them to a single scalar value using parallel processing when enabled.
     * 
     * @tparam BinaryOp The type of the binary operation to apply
     * @param binary_op The binary operation to apply (e.g., std::plus, std::multiplies)
     * @return A scalar value containing the fully reduced result
     * @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
     *       Falls back gracefully to sequential execution if parallel libraries unavailable.
     */
    template<typename BinaryOp>
    T parallel_full_reduce(BinaryOp binary_op) const {
        const size_t total_elements = data_.size();
        return moire_dispatch::reduce_range(data_.begin(), data_.end(), T{}, binary_op, total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD);
    }

    /**
     * @brief Computes the sum of all elements in the multivector
     * Automatically chooses parallel vs sequential based on workload size
     * @return The sum of all elements
     */
    T full_sum() const {
        const size_t total_elements = total_size();
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            return parallel_full_sum();
        }
        // Use sequential for small workloads or when parallelization disabled
        return full_sum(std::execution::seq);
    }

    /**
     * @brief Computes the sum of all elements in the multivector
     * 
     * @tparam ExecutionPolicy The execution policy for parallelization
     * @param policy The execution policy
     * @return The sum of all elements
     */
    template<typename ExecutionPolicy>
    T full_sum(ExecutionPolicy policy) const {
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

    // Smart sum method that automatically chooses parallel vs sequential based on workload size
    MultiVector<T, N-1> sum() const {
        if constexpr (N == 3) {
            const size_t num_slices = dimensions_[0] * dimensions_[1];
            if (should_parallelize(num_slices, MOIRE_PARALLEL_SLICE_THRESHOLD)) {
                return parallel_sum();
            }
        } else if constexpr (N == 2) {
            if (should_parallelize(dimensions_[0], MOIRE_PARALLEL_SLICE_THRESHOLD)) {
                return parallel_sum();
            }
        }
        // Use sequential for small workloads or when parallelization disabled
        return sum(std::execution::seq);
    }

    // Optimized parallel sum method (explicit parallelization)
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
    /// Automatically chooses parallel vs sequential based on workload size
    /// @return A MultiVector<T, N-1> containing the log-sum-exp values
    MultiVector<T, N-1> logsumexp() const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::logsumexp");
#endif
        if constexpr (N == 3) {
            const size_t num_slices = dimensions_[0] * dimensions_[1];
            if (should_parallelize(num_slices, MOIRE_PARALLEL_SLICE_THRESHOLD)) {
                return parallel_logsumexp();
            }
        } else if constexpr (N == 2) {
            if (should_parallelize(dimensions_[0], MOIRE_PARALLEL_SLICE_THRESHOLD)) {
                return parallel_logsumexp();
            }
        }
        // Use sequential for small workloads or when parallelization disabled
        return logsumexp(std::execution::seq);
    }

    /// Compute the log-sum-exp reduction over the innermost dimension
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A MultiVector<T, N-1> containing the log-sum-exp values
    template<typename ExecutionPolicy>
    MultiVector<T, N-1> logsumexp(ExecutionPolicy policy) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::logsumexp(policy)");
#endif
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
                    
                    // Compute exp(x - max) and sum with SIMD vectorization hints
                    T sum = T{};
                    // Use raw pointers with restrict to help compiler optimize
                    const T* __restrict data_ptr = &*begin;
                    const size_t len = static_cast<size_t>(end - begin);
                    
                    // Compiler-specific SIMD vectorization hints
                    #if defined(__GNUC__) || defined(__clang__)
                    #pragma GCC ivdep
                    #elif defined(_MSC_VER)
                    #pragma loop(ivdep)
                    #endif
                    for (size_t k = 0; k < len; ++k) {
                        sum += std::exp(data_ptr[k] - max_val);
                    }
                    
                    // Compute log(sum) + max
                    result.unchecked_at({i, j}) = std::log(sum) + max_val;
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
                
                // Compute exp(x - max) and sum with SIMD vectorization hints
                T sum = T{};
                // Use raw pointers with restrict to help compiler optimize
                const T* __restrict data_ptr = &*begin;
                const size_t len = static_cast<size_t>(end - begin);
                
                // Compiler-specific SIMD vectorization hints
                #if defined(__GNUC__) || defined(__clang__)
                #pragma GCC ivdep
                #elif defined(_MSC_VER)
                #pragma loop(ivdep)
                #endif
                for (size_t k = 0; k < len; ++k) {
                    sum += std::exp(data_ptr[k] - max_val);
                }
                
                // Compute log(sum) + max
                result.unchecked_at({i}) = std::log(sum) + max_val;
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
                    
                    // Compute exp(x - max) and sum with SIMD vectorization hints
                    T sum = T{};
                    // Use raw pointers with restrict to help compiler optimize
                    const T* __restrict data_ptr = &*begin;
                    const size_t len = static_cast<size_t>(end - begin);
                    
                    // Compiler-specific SIMD vectorization hints
                    #if defined(__GNUC__) || defined(__clang__)
                    #pragma GCC ivdep
                    #elif defined(_MSC_VER)
                    #pragma loop(ivdep)
                    #endif
                    for (size_t k = 0; k < len; ++k) {
                        sum += std::exp(data_ptr[k] - max_val);
                    }
                    
                    // Compute log(sum) + max
                    result.unchecked_at(indices) = std::log(sum) + max_val;
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

    /// Compute the log-sum-exp reduction over the innermost dimension with parallel execution when enabled
    /// @return A MultiVector<T, N-1> containing the log-sum-exp values
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    MultiVector<T, N-1> parallel_logsumexp() const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::parallel_logsumexp");
#endif
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
            
            auto logsumexp_slice_fn = [&](size_t i, size_t j) {
                const size_t start_idx = i * stride0 + j * stride1;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim2;
                result.unchecked_at({i, j}) = moire_kernels::logsumexp_slice(begin, end);
            };
            moire_dispatch::for_each_slice_2d(dim0, dim1, MOIRE_PARALLEL_SLICE_THRESHOLD, logsumexp_slice_fn);
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            auto logsumexp_slice_fn = [&](size_t i) {
                const size_t start_idx = i * stride0;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                result.unchecked_at({i}) = moire_kernels::logsumexp_slice(begin, end);
            };
            moire_dispatch::for_each_slice_1d(dim0, MOIRE_PARALLEL_SLICE_THRESHOLD, logsumexp_slice_fn);
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
                    result.unchecked_at(indices) = std::log(sum) + max_val;
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

    /// Apply a transformation to each inner slice of the MultiVector with parallel execution when enabled
    /// @param op The binary operation to apply to each element in the inner slice
    /// @param values A span of values to apply to each inner slice
    /// @return A new MultiVector with the transformed values
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
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
            
            auto transform_slice_2d = [&](size_t i, size_t j) {
                const size_t start_idx = i * stride0 + j * stride1;
                const size_t slice_idx = i * dim1 + j;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim2;
                auto result_begin = result.data_.begin() + start_idx;
                moire_kernels::transform_slice_binary(begin, end, result_begin, values[slice_idx], op);
            };
            moire_dispatch::for_each_slice_2d(dim0, dim1, MOIRE_PARALLEL_SLICE_THRESHOLD, transform_slice_2d);
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            auto transform_slice_1d = [&](size_t i) {
                const size_t start_idx = i * stride0;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                auto result_begin = result.data_.begin() + start_idx;
                moire_kernels::transform_slice_binary(begin, end, result_begin, values[i], op);
            };
            moire_dispatch::for_each_slice_1d(dim0, MOIRE_PARALLEL_SLICE_THRESHOLD, transform_slice_1d);
        }
        else {
            // Fallback to recursive approach for other dimensions
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim, size_t& slice_idx) -> void {
                if (dim == N-1) {
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    std::transform(std::execution::seq, begin, end, result_begin, [&](const T& x) { return op(x, values[slice_idx]); });
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

    /// Apply a transformation to each element within each inner slice with parallel execution when enabled
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element position across all inner slices
    /// @return A new MultiVector with the transformed values
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
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
            
            // Check if we should parallelize (opt-in flag + threshold)
            const size_t total_elements = dim0 * dim1 * dim2;
            if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
                moire_parallel::parallel_for_2d(0, dim0, 0, dim1, [&](size_t i, size_t j) {
                    const size_t start_idx = i * stride0 + j * stride1;
                    for (size_t k = 0; k < dim2; ++k) {
                        result.data_[start_idx + k] = op(data_[start_idx + k], values[k]);
                    }
                });
            } else {
                // Sequential execution when threshold not met or parallelization disabled
                for (size_t i = 0; i < dim0; ++i) {
                    for (size_t j = 0; j < dim1; ++j) {
                        const size_t start_idx = i * stride0 + j * stride1;
                        for (size_t k = 0; k < dim2; ++k) {
                            result.data_[start_idx + k] = op(data_[start_idx + k], values[k]);
                        }
                    }
                }
            }
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
            // Check if we should parallelize (opt-in flag + threshold)
            const size_t total_elements = dim0 * dim1;
            if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
                moire_parallel::parallel_for(0, dim0, [&](size_t i) {
                    const size_t start_idx = i * stride0;
                    for (size_t j = 0; j < dim1; ++j) {
                        result.data_[start_idx + j] = op(data_[start_idx + j], values[j]);
                    }
                });
            } else {
                // Sequential execution when threshold not met or parallelization disabled
                for (size_t i = 0; i < dim0; ++i) {
                    const size_t start_idx = i * stride0;
                    for (size_t j = 0; j < dim1; ++j) {
                        result.data_[start_idx + j] = op(data_[start_idx + j], values[j]);
                    }
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

    // Convenience methods for element-level operations
    /// Add values element-wise, automatically choosing parallel vs sequential based on workload size
    MultiVector<T, N> element_add(const std::span<T const> values) const {
        const size_t total_elements = total_size();
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            return parallel_element_add(values);
        }
        // Use sequential for small workloads or when parallelization disabled
        return element_add(values, std::execution::seq);
    }

    template<typename ExecutionPolicy>
    MultiVector<T, N> element_add(const std::span<T const> values, ExecutionPolicy policy) const {
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
        // Operators are sequential by default to avoid overhead
        // Use parallel_add() if parallelization is needed
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::plus<T>());
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
        // Operators are sequential by default to avoid overhead
        // Use parallel_subtract() if parallelization is needed
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::minus<T>());
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
        // Operators are sequential by default to avoid overhead
        // Use parallel_multiply() if parallelization is needed
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::multiplies<T>());
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
        // Operators are sequential by default to avoid overhead
        // Use parallel_divide() if parallelization is needed
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      result.data_.begin(), std::divides<T>());
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
        // Operators are sequential by default to avoid overhead
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      data_.begin(), std::plus<T>());
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
        // Operators are sequential by default to avoid overhead
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      data_.begin(), std::minus<T>());
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
        // Operators are sequential by default to avoid overhead
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      data_.begin(), std::multiplies<T>());
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
        // Operators are sequential by default to avoid overhead
        std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                      data_.begin(), std::divides<T>());
        return *this;
    }

    /// Parallel elementwise addition of two multivectors with parallel execution when enabled
    /// @param other The other multivector to add
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    MultiVector<T, N> parallel_add(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
        const size_t total_elements = data_.size();
        
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(),
                result.data_.begin(), std::plus<T>());
        } else {
            // Sequential execution when threshold not met or parallelization disabled
            std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                          result.data_.begin(), std::plus<T>());
        }
        return result;
    }

    /// Parallel elementwise subtraction of two multivectors with parallel execution when enabled
    /// @param other The other multivector to subtract
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    MultiVector<T, N> parallel_subtract(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
        const size_t total_elements = data_.size();
        
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(),
                result.data_.begin(), std::minus<T>());
        } else {
            std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                          result.data_.begin(), std::minus<T>());
        }
        return result;
    }

    /// Parallel elementwise multiplication of two multivectors with parallel execution when enabled
    /// @param other The other multivector to multiply with
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    MultiVector<T, N> parallel_multiply(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
        const size_t total_elements = data_.size();
        
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(),
                result.data_.begin(), std::multiplies<T>());
        } else {
            std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                          result.data_.begin(), std::multiplies<T>());
        }
        return result;
    }

    /// Parallel elementwise division of two multivectors with parallel execution when enabled
    /// @param other The other multivector to divide by
    /// @return A new multivector with the result
    /// @throws std::invalid_argument if dimensions don't match
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    MultiVector<T, N> parallel_divide(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        
        MultiVector<T, N> result(dimensions_);
        const size_t total_elements = data_.size();
        
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(),
                result.data_.begin(), std::divides<T>());
        } else {
            std::transform(data_.begin(), data_.end(), other.data_.begin(), 
                          result.data_.begin(), std::divides<T>());
        }
        return result;
    }

    /// Apply softmax transformation to each inner slice of the MultiVector
    /// Automatically chooses parallel vs sequential based on workload size
    /// @param is_log_values Whether the input values are in log space
    /// @return A new MultiVector containing the softmax transformed values
    MultiVector<T, N> softmax(bool is_log_values = true) const {
        if constexpr (N == 3) {
            const size_t num_slices = dimensions_[0] * dimensions_[1];
            if (should_parallelize(num_slices, MOIRE_PARALLEL_SLICE_THRESHOLD)) {
                return parallel_softmax(is_log_values);
            }
        } else if constexpr (N == 2) {
            if (should_parallelize(dimensions_[0], MOIRE_PARALLEL_SLICE_THRESHOLD)) {
                return parallel_softmax(is_log_values);
            }
        }
        // Use sequential for small workloads or when parallelization disabled
        return softmax(is_log_values, std::execution::seq);
    }

    /// Apply softmax transformation to each inner slice of the MultiVector
    /// @param is_log_values Whether the input values are in log space
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A new MultiVector containing the softmax transformed values
    template<typename ExecutionPolicy>
    MultiVector<T, N> softmax(bool is_log_values, ExecutionPolicy policy) const {
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

    /// Apply softmax transformation to each inner slice with parallel execution when enabled
    /// @param is_log_values Whether the input values are in log space
    /// @return A new MultiVector containing the softmax transformed values
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
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
            
            // Check if we should parallelize (opt-in flag + threshold)
            const size_t num_slices = dim0 * dim1;
            if (should_parallelize(num_slices, MOIRE_PARALLEL_SLICE_THRESHOLD)) {
                moire_parallel::parallel_for_2d(0, dim0, 0, dim1, [&](size_t i, size_t j) {
                    const size_t start_idx = i * stride0 + j * stride1;
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    auto result_begin = result.data_.begin() + start_idx;
                    if (is_log_values) {
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        T log_sum = std::log(sum) + max_val;
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = *it - log_sum;
                        }
                    } else {
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = std::exp(*it - max_val) / sum;
                        }
                    }
                });
            } else {
                // Sequential execution when threshold not met or parallelization disabled
                for (size_t i = 0; i < dim0; ++i) {
                    for (size_t j = 0; j < dim1; ++j) {
                        const size_t start_idx = i * stride0 + j * stride1;
                        auto begin = data_.begin() + start_idx;
                        auto end = begin + dim2;
                        auto result_begin = result.data_.begin() + start_idx;
                        
                        if (is_log_values) {
                            T max_val = *std::max_element(std::execution::seq, begin, end);
                            T sum = T{};
                            for (auto it = begin; it != end; ++it) {
                                sum += std::exp(*it - max_val);
                            }
                            T log_sum = std::log(sum) + max_val;
                            for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                                *rit = *it - log_sum;
                            }
                        } else {
                            T max_val = *std::max_element(std::execution::seq, begin, end);
                            T sum = T{};
                            for (auto it = begin; it != end; ++it) {
                                sum += std::exp(*it - max_val);
                            }
                            for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                                *rit = std::exp(*it - max_val) / sum;
                            }
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
            
            // Check if we should parallelize (opt-in flag + threshold)
            if (should_parallelize(dim0, MOIRE_PARALLEL_SLICE_THRESHOLD)) {
                moire_parallel::parallel_for(0, dim0, [&](size_t i) {
                    const size_t start_idx = i * stride0;
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim1;
                    auto result_begin = result.data_.begin() + start_idx;
                    if (is_log_values) {
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        T log_sum = std::log(sum) + max_val;
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = *it - log_sum;
                        }
                    } else {
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = std::exp(*it - max_val) / sum;
                        }
                    }
                });
            } else {
                // Sequential execution when threshold not met or parallelization disabled
                for (size_t i = 0; i < dim0; ++i) {
                    const size_t start_idx = i * stride0;
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim1;
                    auto result_begin = result.data_.begin() + start_idx;
                    
                    if (is_log_values) {
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        T log_sum = std::log(sum) + max_val;
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = *it - log_sum;
                        }
                    } else {
                        T max_val = *std::max_element(std::execution::seq, begin, end);
                        T sum = T{};
                        for (auto it = begin; it != end; ++it) {
                            sum += std::exp(*it - max_val);
                        }
                        for (auto it = begin, rit = result_begin; it != end; ++it, ++rit) {
                            *rit = std::exp(*it - max_val) / sum;
                        }
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
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::transform");
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

    /// Apply a unary transformation to each element with parallel execution when enabled
    /// @param op The unary operation to apply to each element
    /// @return A new MultiVector with the transformed values
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    template<typename UnaryOp>
    MultiVector<T, N> parallel_transform(UnaryOp op) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::parallel_transform");
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
            
            // Check if we should parallelize (opt-in flag + threshold)
            const size_t total_elements = dim0 * dim1 * dim2;
            if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
                moire_parallel::parallel_for_2d(0, dim0, 0, dim1, [&](size_t i, size_t j) {
                    const size_t start_idx = i * stride0 + j * stride1;
                    auto begin = data_.begin() + start_idx;
                    auto end = begin + dim2;
                    auto result_begin = result.data_.begin() + start_idx;
                    std::transform(std::execution::seq, begin, end, result_begin, op);
                });
            } else {
                // Sequential execution when threshold not met or parallelization disabled
                std::transform(data_.begin(), data_.end(), result.data_.begin(), op);
            }
        }
        else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            
            // Pre-calculate stride for faster access
            const size_t stride0 = strides_[0];
            
            // Check if we should parallelize (opt-in flag + threshold)
            const size_t total_elements = dim0 * dim1;
            if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
                moire_parallel::transform(data_.begin(), data_.end(), result.data_.begin(), op);
            } else {
                // Sequential execution when threshold not met or parallelization disabled
                std::transform(data_.begin(), data_.end(), result.data_.begin(), op);
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

    /// Compute natural logarithm of each element
    /// Automatically chooses parallel vs sequential based on workload size
    /// @return A MultiVector with log-transformed values
    MultiVector<T, N> log() const {
        const size_t total_elements = total_size();
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            return parallel_log();
        }
        // Use sequential for small workloads or when parallelization disabled
        return log(std::execution::seq);
    }

    template<typename ExecutionPolicy>
    MultiVector<T, N> log(ExecutionPolicy policy) const {
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

    std::vector<T> data_;
    std::array<size_t, N> dimensions_;
    std::array<size_t, N> strides_;

    /// Unchecked element access for internal hot paths. Caller must ensure indices are valid.
    inline T& unchecked_at(const std::array<size_t, N>& indices) {
        return data_[calculate_index(indices)];
    }
    inline const T& unchecked_at(const std::array<size_t, N>& indices) const {
        return data_[calculate_index(indices)];
    }

    // Helper function to calculate the 1D index from multi-dimensional indices
    /// @param indices The multi-dimensional indices.
    /// @return The calculated 1D index.
    /// @throws std::out_of_range if any index is out of bounds.
    inline size_t calculate_index(const std::array<size_t, N>& indices) const {
#ifndef NDEBUG
        for (size_t i = 0; i < N; ++i) {
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

#include "multivector_1d.h"
#include "ragged_multivector.h"

