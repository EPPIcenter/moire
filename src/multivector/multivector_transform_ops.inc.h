// Slice-level transforms (logsumexp, binary slice ops, softmax, unary) for MultiVector<T,N>.
// Included inside the MultiVector class body from multivector.h.

    MultiVector<T, N-1> logsumexp() const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::logsumexp");
#endif
        std::array<size_t, N-1> reduced_dims;
        std::copy(dimensions_.begin(), dimensions_.end() - 1, reduced_dims.begin());
        MultiVector<T, N-1> result(reduced_dims);

        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            moire_parallel::for_each_slice_2d(dim0, dim1, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i, size_t j) {
                const size_t start_idx = i * stride0 + j * stride1;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim2;
                result.unchecked_at({i, j}) = moire_kernels::logsumexp_slice(begin, end);
            });
        } else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            moire_parallel::for_each_slice_1d(dim0, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i) {
                const size_t start_idx = i * stride0;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                result.unchecked_at({i}) = moire_kernels::logsumexp_slice(begin, end);
            });
        } else {
            auto logsumexp_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    const auto [begin, end] = inner_iterators(indices);
                    result.unchecked_at(indices) = moire_kernels::logsumexp_slice(begin, end);
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

    template<typename BinaryOp>
    MultiVector<T, N> transform(BinaryOp op, const std::span<T const> values) const {
        MultiVector<T, N> result(dimensions_);
        size_t num_inner_slices = 1;
        for (size_t i = 0; i < N - 1; ++i) num_inner_slices *= dimensions_[i];
        if (values.size() != num_inner_slices) {
            throw std::invalid_argument("Values span size does not match the number of inner slices");
        }

        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            moire_parallel::for_each_slice_2d(dim0, dim1, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i, size_t j) {
                const size_t start_idx = i * stride0 + j * stride1;
                const size_t slice_idx = i * dim1 + j;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim2;
                auto result_begin = result.data_.begin() + start_idx;
                moire_kernels::transform_slice_binary(begin, end, result_begin, values[slice_idx], op);
            });
        } else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            moire_parallel::for_each_slice_1d(dim0, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i) {
                const size_t start_idx = i * stride0;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                auto result_begin = result.data_.begin() + start_idx;
                moire_kernels::transform_slice_binary(begin, end, result_begin, values[i], op);
            });
        } else {
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim, size_t& slice_idx) -> void {
                if (dim == N-1) {
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    moire_kernels::transform_slice_binary(begin, end, result_begin, values[slice_idx], op);
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

    MultiVector<T, N> add(const std::span<T const> values) const { return transform(std::plus<T>(), values); }
    MultiVector<T, N> subtract(const std::span<T const> values) const { return transform(std::minus<T>(), values); }
    MultiVector<T, N> multiply(const std::span<T const> values) const { return transform(std::multiplies<T>(), values); }
    MultiVector<T, N> divide(const std::span<T const> values) const { return transform(std::divides<T>(), values); }

    MultiVector<T, N> softmax(bool is_log_values = true) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::softmax");
#endif
        MultiVector<T, N> result(dimensions_);

        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            moire_parallel::for_each_slice_2d(dim0, dim1, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i, size_t j) {
                const size_t start_idx = i * stride0 + j * stride1;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim2;
                auto result_begin = result.data_.begin() + start_idx;
                moire_kernels::softmax_slice(begin, end, result_begin, is_log_values);
            });
        } else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            moire_parallel::for_each_slice_1d(dim0, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i) {
                const size_t start_idx = i * stride0;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                auto result_begin = result.data_.begin() + start_idx;
                moire_kernels::softmax_slice(begin, end, result_begin, is_log_values);
            });
        } else {
            auto softmax_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    moire_kernels::softmax_slice(begin, end, result_begin, is_log_values);
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

    template<typename UnaryOp>
    MultiVector<T, N> transform(UnaryOp op) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::transform");
#endif
        MultiVector<T, N> result(dimensions_);
        if constexpr (N == 3) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t dim2 = dimensions_[2];
            const size_t stride0 = strides_[0];
            const size_t stride1 = strides_[1];
            const size_t total_elements = dim0 * dim1 * dim2;
            if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)
                && dim2 >= MOIRE_PARALLEL_SLICE_UNARY_INNER_THRESHOLD) {
                moire_parallel::parallel_for_2d(0, dim0, 0, dim1, [&](size_t i, size_t j) {
                    const size_t start_idx = i * stride0 + j * stride1;
                    moire_kernels::transform_slice_unary(
                        data_.begin() + start_idx, data_.begin() + start_idx + dim2,
                        result.data_.begin() + start_idx, op);
                });
            } else if (should_parallelize(total_elements, MOIRE_PARALLEL_ELEMENT_THRESHOLD)) {
                moire_parallel::transform(data_.begin(), data_.end(), result.data_.begin(), op);
            } else {
                moire_kernels::transform_slice_unary(data_.begin(), data_.end(), result.data_.begin(), op);
            }
        } else if constexpr (N == 2) {
            moire_parallel::transform(data_.begin(), data_.end(), result.data_.begin(), op);
        } else {
            auto transform_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    const auto [begin, end] = inner_iterators(indices);
                    const auto [result_begin, _] = result.inner_iterators(indices);
                    moire_kernels::transform_slice_unary(begin, end, result_begin, op);
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

    MultiVector<T, N> negate() const { return transform([](const T& x) { return -x; }); }
    MultiVector<T, N> abs() const { return transform([](const T& x) { return std::abs(x); }); }
    MultiVector<T, N> exp() const { return transform([](const T& x) { return std::exp(x); }); }
    MultiVector<T, N> log() const { return transform([](const T& x) { return std::log(x); }); }
    MultiVector<T, N> sqrt() const { return transform([](const T& x) { return std::sqrt(x); }); }
