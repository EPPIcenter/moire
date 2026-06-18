#pragma once

// MultiVector<T,1> specialization; include after primary MultiVector template in multivector.h

template <typename T>
class MultiVector<T, 1> {
public:
    MultiVector(const std::array<size_t, 1>& dimensions) : dimensions_(dimensions) {
        strides_[0] = 1;
        data_.resize(dimensions_[0], T{});
    }

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

    std::span<T> as_span() { return std::span<T>(data_); }
    std::span<const T> as_span() const { return std::span<const T>(data_); }
    operator std::span<T>() { return as_span(); }
    operator std::span<const T>() const { return as_span(); }

    void inner_fill(const T& value) { std::fill(data_.begin(), data_.end(), value); }

    void inner_fill(const std::span<T const> values) {
        if (values.size() > dimensions_[0]) {
            throw std::invalid_argument("Incorrect number of values for 1D inner_fill");
        }
        std::copy(values.begin(), values.end(), data_.begin());
    }

    T& at(const std::array<size_t, 1>& indices) {
        if (indices[0] >= dimensions_[0]) throw std::out_of_range("MultiVector::at index out of bounds");
        return data_.at(indices[0]);
    }

    const T& at(const std::array<size_t, 1>& indices) const {
        if (indices[0] >= dimensions_[0]) throw std::out_of_range("MultiVector::at index out of bounds");
        return data_.at(indices[0]);
    }

    size_t size(size_t dimension) const {
        if (dimension >= 1u) throw std::out_of_range("Dimension out of range");
        return dimensions_[dimension];
    }

    size_t total_size() const { return data_.size(); }

    typename std::vector<T>::iterator inner_begin() { return data_.begin(); }
    typename std::vector<T>::iterator inner_end() { return data_.end(); }
    typename std::vector<T>::const_iterator inner_begin() const { return data_.cbegin(); }
    typename std::vector<T>::const_iterator inner_end() const { return data_.cend(); }
    std::pair<typename std::vector<T>::const_iterator, typename std::vector<T>::const_iterator> inner_iterators() const {
        return {data_.cbegin(), data_.cend()};
    }

    void clear() { data_.clear(); }
    std::array<size_t, 1> dimensions() const { return dimensions_; }
    const std::array<size_t, 1>& strides() const { return strides_; }
    const std::vector<T>& data() const { return data_; }

    template<typename BinaryOp>
    auto reduce(BinaryOp binary_op, const auto& init) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector<T,1>::reduce");
#endif
        return moire_parallel::reduce_range(
            data_.begin(), data_.end(), init, binary_op, data_.size(), MOIRE_PARALLEL_ELEMENT_THRESHOLD);
    }

    T full_sum() const { return sum(); }
    T full_product() const { return product(); }
    T full_min() const { return min(); }
    T full_max() const { return max(); }

    T sum() const { return reduce(std::plus<T>(), T{}); }
    T product() const { return reduce(std::multiplies<T>(), T{1}); }
    T max() const { return reduce([](const T& a, const T& b) { return std::max(a, b); }, std::numeric_limits<T>::lowest()); }
    T min() const { return reduce([](const T& a, const T& b) { return std::min(a, b); }, std::numeric_limits<T>::max()); }

    T logsumexp() const { return moire_kernels::logsumexp_slice(data_.begin(), data_.end()); }

    MultiVector<T, 1> element_add(const std::span<T const> values) const { return transform(std::plus<T>(), values); }
    MultiVector<T, 1> element_subtract(const std::span<T const> values) const { return transform(std::minus<T>(), values); }
    MultiVector<T, 1> element_multiply(const std::span<T const> values) const { return transform(std::multiplies<T>(), values); }
    MultiVector<T, 1> element_divide(const std::span<T const> values) const { return transform(std::divides<T>(), values); }

    template<typename UnaryOp>
    MultiVector<T, 1> transform(UnaryOp op) const {
        MultiVector<T, 1> result(dimensions_);
        moire_parallel::transform(data_.begin(), data_.end(), result.data_.begin(), op);
        return result;
    }

    MultiVector<T, 1> negate() const { return transform([](const T& x) { return -x; }); }
    MultiVector<T, 1> abs() const { return transform([](const T& x) { return std::abs(x); }); }
    MultiVector<T, 1> exp() const { return transform([](const T& x) { return std::exp(x); }); }
    MultiVector<T, 1> log() const { return transform([](const T& x) { return std::log(x); }); }
    MultiVector<T, 1> sqrt() const { return transform([](const T& x) { return std::sqrt(x); }); }

    template<typename BinaryOp>
    MultiVector<T, 1> transform(BinaryOp op, const std::span<T const> values) const {
        if (values.size() != dimensions_[0]) {
            throw std::invalid_argument("Values span size does not match the vector size");
        }
        MultiVector<T, 1> result(dimensions_);
        moire_parallel::transform(data_.begin(), data_.end(), values.begin(), result.data_.begin(), op);
        return result;
    }

    MultiVector<T, 1> add(const std::span<T const> values) const { return transform(std::plus<T>(), values); }
    MultiVector<T, 1> subtract(const std::span<T const> values) const { return transform(std::minus<T>(), values); }
    MultiVector<T, 1> multiply(const std::span<T const> values) const { return transform(std::multiplies<T>(), values); }
    MultiVector<T, 1> divide(const std::span<T const> values) const { return transform(std::divides<T>(), values); }

    std::vector<T> data_;
    std::array<size_t, 1> dimensions_;
    std::array<size_t, 1> strides_;

    inline T& unchecked_at(const std::array<size_t, 1>& indices) { return data_[calculate_index(indices)]; }
    inline const T& unchecked_at(const std::array<size_t, 1>& indices) const { return data_[calculate_index(indices)]; }

    inline size_t calculate_index(const std::array<size_t, 1>& indices) const { return indices[0]; }

    std::pair<size_t, size_t> calculate_start_end_indices(const std::array<size_t, 0>&) const {
        return {0, dimensions_[0]};
    }
};
