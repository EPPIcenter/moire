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

    /// Invoke func(outer_indices, begin, end) for each innermost slice.
    template<typename Func>
    void for_each_inner_slice(Func&& func) {
        if constexpr (N == 1) {
            std::array<size_t, 0> outer{};
            func(outer, data_.begin(), data_.end());
        } else if constexpr (N == 2) {
            for (size_t i = 0; i < dimensions_[0]; ++i) {
                std::array<size_t, 1> outer{i};
                auto [begin, end] = inner_iterators(outer);
                func(outer, begin, end);
            }
        } else if constexpr (N == 3) {
            for (size_t i = 0; i < dimensions_[0]; ++i) {
                for (size_t j = 0; j < dimensions_[1]; ++j) {
                    std::array<size_t, 2> outer{i, j};
                    auto [begin, end] = inner_iterators(outer);
                    func(outer, begin, end);
                }
            }
        } else {
            auto visit = [&](auto& self, std::array<size_t, N - 1>& outer, size_t dim) -> void {
                if (dim == N - 1) {
                    auto [begin, end] = inner_iterators(outer);
                    func(outer, begin, end);
                    return;
                }
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    outer[dim] = i;
                    self(self, outer, dim + 1);
                }
            };
            std::array<size_t, N - 1> outer{};
            visit(visit, outer, 0);
        }
    }

    template<typename Func>
    void for_each_inner_slice(Func&& func) const {
        if constexpr (N == 1) {
            std::array<size_t, 0> outer{};
            func(outer, data_.cbegin(), data_.cend());
        } else if constexpr (N == 2) {
            for (size_t i = 0; i < dimensions_[0]; ++i) {
                std::array<size_t, 1> outer{i};
                auto [begin, end] = inner_iterators(outer);
                func(outer, begin, end);
            }
        } else if constexpr (N == 3) {
            for (size_t i = 0; i < dimensions_[0]; ++i) {
                for (size_t j = 0; j < dimensions_[1]; ++j) {
                    std::array<size_t, 2> outer{i, j};
                    auto [begin, end] = inner_iterators(outer);
                    func(outer, begin, end);
                }
            }
        } else {
            auto visit = [&](auto& self, std::array<size_t, N - 1>& outer, size_t dim) -> void {
                if (dim == N - 1) {
                    auto [begin, end] = inner_iterators(outer);
                    func(outer, begin, end);
                    return;
                }
                for (size_t i = 0; i < dimensions_[dim]; ++i) {
                    outer[dim] = i;
                    self(self, outer, dim + 1);
                }
            };
            std::array<size_t, N - 1> outer{};
            visit(visit, outer, 0);
        }
    }


#include "multivector_reduce_ops.inc.h"
#include "multivector_transform_ops.inc.h"
#include "multivector_element_ops.inc.h"

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

