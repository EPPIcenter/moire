// Element-wise and MultiVector-MultiVector operations for MultiVector<T,N>.
// Included inside the MultiVector class body from multivector.h.

    template<typename BinaryOp>
    MultiVector<T, N> element_transform(BinaryOp op, const std::span<T const> values) const {
        if (values.size() != dimensions_.back()) {
            throw std::invalid_argument("Values span size does not match the inner dimension size");
        }

        MultiVector<T, N> result(dimensions_);
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            moire_parallel::for_each_slice_2d(dim0, dim1, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i, size_t j) {
                const size_t start_idx = i * stride0 + j * stride1;
                for (size_t k = 0; k < dim2; ++k) {
                    result.data_[start_idx + k] = op(data_[start_idx + k], values[k]);
                }
            });
        } else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            moire_parallel::parallel_for(0, dim0, [&](size_t i) {
                const size_t start_idx = i * stride0;
                for (size_t j = 0; j < dim1; ++j) {
                    result.data_[start_idx + j] = op(data_[start_idx + j], values[j]);
                }
            });
        } else {
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
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

    MultiVector<T, N> element_add(const std::span<T const> values) const {
        return element_transform(std::plus<T>(), values);
    }
    MultiVector<T, N> element_subtract(const std::span<T const> values) const {
        return element_transform(std::minus<T>(), values);
    }
    MultiVector<T, N> element_multiply(const std::span<T const> values) const {
        return element_transform(std::multiplies<T>(), values);
    }
    MultiVector<T, N> element_divide(const std::span<T const> values) const {
        return element_transform(std::divides<T>(), values);
    }

    bool has_same_dimensions(const MultiVector<T, N>& other) const {
        return dimensions_ == other.dimensions_;
    }

    MultiVector<T, N> operator+(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        MultiVector<T, N> result(dimensions_);
        moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(), result.data_.begin(), std::plus<T>());
        return result;
    }

    MultiVector<T, N> operator-(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        MultiVector<T, N> result(dimensions_);
        moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(), result.data_.begin(), std::minus<T>());
        return result;
    }

    MultiVector<T, N> operator*(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        MultiVector<T, N> result(dimensions_);
        moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(), result.data_.begin(), std::multiplies<T>());
        return result;
    }

    MultiVector<T, N> operator/(const MultiVector<T, N>& other) const {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        MultiVector<T, N> result(dimensions_);
        moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(), result.data_.begin(), std::divides<T>());
        return result;
    }

    MultiVector<T, N>& operator+=(const MultiVector<T, N>& other) {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(), data_.begin(), std::plus<T>());
        return *this;
    }

    MultiVector<T, N>& operator-=(const MultiVector<T, N>& other) {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(), data_.begin(), std::minus<T>());
        return *this;
    }

    MultiVector<T, N>& operator*=(const MultiVector<T, N>& other) {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(), data_.begin(), std::multiplies<T>());
        return *this;
    }

    MultiVector<T, N>& operator/=(const MultiVector<T, N>& other) {
        if (!has_same_dimensions(other)) {
            throw std::invalid_argument("Multivectors must have the same dimensions for elementwise operations");
        }
        moire_parallel::transform(data_.begin(), data_.end(), other.data_.begin(), data_.begin(), std::divides<T>());
        return *this;
    }
