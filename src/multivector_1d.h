#pragma once

// MultiVector<T,1> specialization; include after primary MultiVector template in multivector.h

/// Specialization of MultiVector for N=1 that returns a scalar from reduce functions
template <typename T>
class MultiVector<T, 1> {
public:
    /// Constructor that takes dimensions
    /// @param dimensions The dimensions of the MultiVector.
    MultiVector(const std::array<size_t, 1>& dimensions) : dimensions_(dimensions) {
        strides_[0] = 1;
        data_.resize(dimensions_[0], T{});
    }

    /// Default constructor: creates an empty (valid) 1D MultiVector with total_size() == 0.
    MultiVector() {
        dimensions_ = {0};
        strides_[0] = 1;
        data_.clear();
    }

    void resize(const std::array<size_t, 1>& dimensions) {
        dimensions_ = dimensions;
        strides_[0] = 1;
        data_.resize(dimensions_[0], T{});
    }

    void resize(const std::array<size_t, 1>& dimensions, const T& value) {
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

    /// Access element at a given index (always bounds-checked).
    T& at(const std::array<size_t, 1>& indices) {
        if (indices[0] >= dimensions_[0]) {
            throw std::out_of_range("MultiVector::at index out of bounds");
        }
        return data_.at(indices[0]);
    }

    const T& at(const std::array<size_t, 1>& indices) const {
        if (indices[0] >= dimensions_[0]) {
            throw std::out_of_range("MultiVector::at index out of bounds");
        }
        return data_.at(indices[0]);
    }

    /// Get the size of the dimension
    /// @return The size of the dimension.
    size_t size(size_t dimension) const {
        if (dimension >= 1u) {
            throw std::out_of_range("Dimension out of range");
        }
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

    /// Optimized parallel reduce using parallel execution when enabled and workload exceeds threshold
    /// @param binary_op The binary operation to use for reduction
    /// @param init The initial value for the reduction
    /// @return A scalar value containing the reduced value
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    template<typename BinaryOp>
    auto parallel_reduce(BinaryOp binary_op, const auto& init) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector<T,1>::parallel_reduce");
#endif
        const size_t total_elements = data_.size();
        return moire_dispatch::reduce_range(data_.begin(), data_.end(), init, binary_op, total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD);
    }

    // Convenience methods for common reductions with optional execution policy
    template<typename ExecutionPolicy = std::execution::sequenced_policy>
    T sum(ExecutionPolicy policy = std::execution::seq) const {
        return reduce(std::plus<T>(), T{}, policy);
    }

    // Smart sum method that automatically chooses parallel vs sequential based on workload size
    T sum() const {
        const size_t total_elements = total_size();
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            return parallel_sum();
        }
        // Use sequential for small workloads or when parallelization disabled
        return sum(std::execution::seq);
    }

    // Optimized parallel sum method (explicit parallelization)
    T parallel_sum() const {
        return parallel_reduce(std::plus<T>(), T{});
    }

    T full_sum() const {
        return sum();  // Will use smart sum() that chooses automatically
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
    /// Automatically chooses parallel vs sequential based on workload size
    /// @return A scalar value containing the log-sum-exp value
    T logsumexp() const {
        const size_t total_elements = total_size();
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            return parallel_logsumexp();
        }
        // Use sequential for small workloads or when parallelization disabled
        return logsumexp(std::execution::seq);
    }

    /// Compute the log-sum-exp reduction over the vector
    /// @param execution_policy Optional execution policy for parallel execution
    /// @return A scalar value containing the log-sum-exp value
    template<typename ExecutionPolicy>
    T logsumexp(ExecutionPolicy policy) const {
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

    /// Compute the log-sum-exp reduction over the vector with parallel execution when enabled
    /// @return A scalar value containing the log-sum-exp value
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    T parallel_logsumexp() const {
        const size_t total_elements = data_.size();
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            T max_val = moire_parallel::reduce(data_.begin(), data_.end(), *data_.begin(),
                [](const T& a, const T& b) { return std::max(a, b); });
            T sum = moire_parallel::transform_reduce(data_.begin(), data_.end(), T{}, std::plus<T>(),
                [max_val](const T& x) { return std::exp(x - max_val); });
            return std::log(sum) + max_val;
        }
        return moire_kernels::logsumexp_slice(data_.begin(), data_.end());
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

    /// Apply a transformation to each element of the vector with parallel execution when enabled
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    template<typename BinaryOp>
    void parallel_element_transform(BinaryOp op, const std::span<T const> values) {
#ifndef NDEBUG
        if (values.size() != dimensions_[0]) {
            throw std::invalid_argument("Values span size does not match the vector size");
        }
#endif
        const size_t total_elements = dimensions_[0];
        
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            moire_parallel::transform(data_.begin(), data_.end(), values.begin(), data_.begin(), op);
        } else {
            // Sequential execution when threshold not met or parallelization disabled
            std::transform(std::execution::seq, data_.begin(), data_.end(), values.begin(), data_.begin(), op);
        }
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

    inline T& unchecked_at(const std::array<size_t, 1>& indices) {
        return data_[calculate_index(indices)];
    }
    inline const T& unchecked_at(const std::array<size_t, 1>& indices) const {
        return data_[calculate_index(indices)];
    }

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

    /// Apply a unary transformation to each element with parallel execution when enabled
    /// @param op The unary operation to apply to each element
    /// @return A new MultiVector with the transformed values
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    template<typename UnaryOp>
    MultiVector<T, 1> parallel_transform(UnaryOp op) const {
        MultiVector<T, 1> result(dimensions_);
        const size_t total_elements = dimensions_[0];
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            moire_parallel::transform(data_.begin(), data_.end(), result.data_.begin(), op);
        } else {
            std::transform(std::execution::seq, data_.begin(), data_.end(), result.data_.begin(), op);
        }
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

    /// Compute natural logarithm of each element
    /// Automatically chooses parallel vs sequential based on workload size
    /// @return A MultiVector with log-transformed values
    MultiVector<T, 1> log() const {
        const size_t total_elements = total_size();
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            return parallel_log();
        }
        // Use sequential for small workloads or when parallelization disabled
        return log(std::execution::seq);
    }

    template<typename ExecutionPolicy>
    MultiVector<T, 1> log(ExecutionPolicy policy) const {
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

    /// Apply a transformation to each element with parallel execution when enabled
    /// @param op The binary operation to apply to each element
    /// @param values A span of values to apply to each element
    /// @return A new MultiVector with the transformed values
    /// @note Requires MOIRE_ENABLE_PARALLEL to be defined for parallel execution.
    ///       Falls back gracefully to sequential execution if parallel libraries unavailable.
    template<typename BinaryOp>
    MultiVector<T, 1> parallel_transform(BinaryOp op, const std::span<T const> values) const {
        if (values.size() != dimensions_[0]) {
            throw std::invalid_argument("Values span size does not match the vector size");
        }
        MultiVector<T, 1> result(dimensions_);
        const size_t total_elements = dimensions_[0];
        if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
            moire_parallel::transform(data_.begin(), data_.end(), values.begin(), result.data_.begin(), op);
        } else {
            std::transform(std::execution::seq, data_.begin(), data_.end(), values.begin(), result.data_.begin(), op);
        }
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

