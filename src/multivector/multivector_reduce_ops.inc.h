// Slice-reduction operations for MultiVector<T,N>.
// Included inside the MultiVector class body from multivector.h.

    template<typename BinaryOp>
    auto reduce(BinaryOp binary_op, const auto& init) const {
#ifdef MOIRE_ENABLE_PROFILER_REGISTRY
        ProfileScope _prof("MultiVector::reduce");
#endif
        using ResultType = decltype(binary_op(init, std::declval<T>()));
        std::array<size_t, N-1> reduced_dims;
        std::copy(dimensions_.begin(), dimensions_.end() - 1, reduced_dims.begin());
        MultiVector<ResultType, N-1> result(reduced_dims);

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
                result.unchecked_at({i, j}) = moire_kernels::reduce_slice(begin, end, init, binary_op);
            });
        } else if constexpr (N == 2) {
            const size_t dim0 = dimensions_[0];
            const size_t dim1 = dimensions_[1];
            const size_t stride0 = strides_[0];
            moire_parallel::for_each_slice_1d(dim0, MOIRE_PARALLEL_SLICE_THRESHOLD, [&](size_t i) {
                const size_t start_idx = i * stride0;
                auto begin = data_.begin() + start_idx;
                auto end = begin + dim1;
                result.unchecked_at({i}) = moire_kernels::reduce_slice(begin, end, init, binary_op);
            });
        } else {
            auto reduce_recursive = [&](auto& self, std::array<size_t, N-1>& indices, size_t dim) -> void {
                if (dim == N-1) {
                    const auto [begin, end] = inner_iterators(indices);
                    result.unchecked_at(indices) = moire_kernels::reduce_slice(begin, end, init, binary_op);
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

    template<typename BinaryOp>
    T full_reduce(BinaryOp binary_op) const {
        return moire_parallel::reduce_range(
            data_.begin(), data_.end(), T{}, binary_op, data_.size(), MOIRE_PARALLEL_ELEMENT_THRESHOLD);
    }

    T full_sum() const { return full_reduce(std::plus<T>()); }
    T full_product() const { return full_reduce(std::multiplies<T>()); }
    T full_min() const { return full_reduce([](const T& a, const T& b) { return std::min(a, b); }); }
    T full_max() const { return full_reduce([](const T& a, const T& b) { return std::max(a, b); }); }

    MultiVector<T, N-1> sum() const { return reduce(std::plus<T>(), T{}); }
    MultiVector<T, N-1> product() const { return reduce(std::multiplies<T>(), T{1}); }
    MultiVector<T, N-1> max() const {
        return reduce([](const T& a, const T& b) { return std::max(a, b); }, std::numeric_limits<T>::lowest());
    }
    MultiVector<T, N-1> min() const {
        return reduce([](const T& a, const T& b) { return std::min(a, b); }, std::numeric_limits<T>::max());
    }
