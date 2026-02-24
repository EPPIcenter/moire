#pragma once

#include <vector>
#include <array>
#include <span>
#include <stdexcept>
#include <execution>
#include <numeric>
#include <algorithm>
#include <cmath>

/// RaggedMultiVector class template (N >= 2).
/// Multi-dimensional vector with a ragged innermost dimension per outer slot.
///
/// Canonical contract:
/// - dimensions_.size() == N-1 (outer shape). total_entries = product(dimensions_[0..N-3]).
/// - ragged_dimensions_.size() is the number of ragged slots per block; each block has the same layout.
/// - data_.size() == total_entries * sum(ragged_dimensions_). ragged_offsets_[i] = sum(ragged_dimensions_[0..i-1]); ragged_offsets_[0] == 0.
/// - dimensions() returns [dimensions_[0], ..., dimensions_[N-2], ragged_dimensions_.size()]; the last component is the number of slots per block.
template <typename T, size_t N>
class RaggedMultiVector {
    static_assert(N >= 2, "RaggedMultiVector requires N >= 2");
public:
    /// Constructor that takes dimensions and ragged dimensions.
    /// @param dimensions Outer shape (N-1 elements).
    /// @param ragged_dimensions Length of innermost dimension per slot; one block has ragged_dimensions_.size() slots.
    /// @param default_value The default value for elements.
    RaggedMultiVector(const std::array<size_t, N - 1>& dimensions, const std::span<size_t const> ragged_dimensions, const T& default_value = T{}) : dimensions_(dimensions), ragged_dimensions_(ragged_dimensions.begin(), ragged_dimensions.end()) {
        std::vector<size_t> entry_dimensions(dimensions_.begin(), dimensions_.end() - 1);
        size_t total_entries = std::reduce(std::execution::seq, entry_dimensions.begin(), entry_dimensions.end(), size_t(1), std::multiplies<size_t>());
        {
            size_t data_size = std::reduce(std::execution::seq, ragged_dimensions_.begin(), ragged_dimensions_.end());
            size_t total_size = total_entries * data_size;

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


    /// Default constructor: creates an empty (valid) RaggedMultiVector.
    RaggedMultiVector() {
        dimensions_ = {};
        ragged_dimensions_.clear();
        ragged_offsets_.clear();
        data_.clear();
    }

    void resize(const std::array<size_t, N - 1>& dimensions, const std::span<size_t const> ragged_dimensions, const T& default_value = T{}) {
        dimensions_ = dimensions;
        ragged_dimensions_ = std::vector<size_t>(ragged_dimensions.begin(), ragged_dimensions.end());
        std::vector<size_t> entry_dimensions(dimensions_.begin(), dimensions_.end() - 1);
        size_t total_entries = std::reduce(std::execution::seq, entry_dimensions.begin(), entry_dimensions.end(), size_t(1), std::multiplies<size_t>());
        size_t data_size = std::reduce(std::execution::seq, ragged_dimensions_.begin(), ragged_dimensions_.end());
        size_t total_size = total_entries * data_size;

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
        if (indices.back() >= ragged_dimensions_.size()) {
            throw std::out_of_range("Index out of bounds for ragged_dimensions_");
        }
        if (values.size() > ragged_dimensions_.at(indices.back())) {
            throw std::invalid_argument("(" + std::to_string(N) + "D) Incorrect number of values, expected " + std::to_string(ragged_dimensions_.at(indices.back())) + " but got " + std::to_string(values.size()));
        }
        std::copy(values.begin(), values.end(), inner_begin(indices));
    }

    size_t total_size() const {
        return data_.size();
    }

    /// Access element at a given multi-dimensional index (always bounds-checked).
    T& at(const std::array<size_t, N>& indices) {
        for (size_t i = 0; i < N - 1; ++i) {
            if (indices[i] >= dimensions_[i]) {
                throw std::out_of_range("RaggedMultiVector::at index out of bounds");
            }
        }
        if (indices[N - 2] >= ragged_dimensions_.size() || indices.back() >= ragged_dimensions_.at(indices[N - 2])) {
            throw std::out_of_range("RaggedMultiVector::at index out of bounds");
        }
        return data_.at(calculate_index(indices));
    }

    const T& at(const std::array<size_t, N>& indices) const {
        for (size_t i = 0; i < N - 1; ++i) {
            if (indices[i] >= dimensions_[i]) {
                throw std::out_of_range("RaggedMultiVector::at index out of bounds");
            }
        }
        if (indices[N - 2] >= ragged_dimensions_.size() || indices.back() >= ragged_dimensions_.at(indices[N - 2])) {
            throw std::out_of_range("RaggedMultiVector::at index out of bounds");
        }
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

private:
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
        for (size_t i = 0; i < N - 1; ++i) {
            if (indices[i] >= dimensions_[i]) {
                throw std::out_of_range("Index out of bounds");
            }
        }
        if (indices.back() >= ragged_dimensions_.at(indices[N - 2])) {
            throw std::out_of_range("Index out of bounds");
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
        std::array<size_t, N> full_indices{};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        const size_t start_index = calculate_index(full_indices);
        const size_t end_index = start_index + ragged_dimensions_.at(outer_indices[N - 2]);
        return {start_index, end_index};
    }
};

