#pragma once

#include <vector>
#include <stdexcept>
#include <iostream>
#include <numeric>
#include <array>

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
    typename std::span<T>::iterator inner_begin(const std::array<size_t, N-1>& outer_indices) {
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

    typename std::span<T>::iterator inner_end(const std::array<size_t, N-1>& outer_indices) {
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

    /// Clear the MultiVector
    /// @note This will clear the data contained in the MultiVector.
    void clear() {
        data_.clear();
    }

    std::array<size_t, N> dimensions() const {
        return dimensions_;
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
        size_t index = 0;
        const size_t* stride = strides_.data();
        const size_t* idx = indices.data();
        for (size_t i = N; i-- > 0;) {
            index += idx[i] * stride[i];
        }
        return index;
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
    RaggedMultiVector(const std::array<size_t, N - 1>& dimensions, const std::span<size_t const> ragged_dimensions, const T& default_value = T{}) : dimensions_(dimensions), ragged_dimensions_(ragged_dimensions) {
        // Calculate total size and strides
        size_t total_size = 0;
        for (const auto& rd : ragged_dimensions_) {
            total_size += rd;
        }
        for (size_t i = N - 1; i-- > 0;) {
            strides_[i] = total_size;
            total_size *= dimensions_[i];
        }
        // ragged strides are the cumulative sum of the ragged dimensions and provide the offset necessary to access the ragged dimension.
        ragged_offsets_.resize(ragged_dimensions_.size() + 1, 0);
        for (size_t i = 0; i < ragged_dimensions_.size(); ++i) {
            ragged_offsets_.at(i + 1) = ragged_offsets_.at(i) + ragged_dimensions_.at(i);
        }
        data_.resize(total_size, default_value);
    }

    RaggedMultiVector() {
        dimensions_ = {0};
        ragged_dimensions_ = {};
        ragged_offsets_ = {};
        data_ = {};
    }

    void resize(const std::array<size_t, N - 1>& dimensions, std::span<size_t const> ragged_dimensions, const T& default_value = T{}) {
        dimensions_ = dimensions;
        ragged_dimensions_ = std::vector<size_t>(ragged_dimensions.begin(), ragged_dimensions.end());
        size_t total_size = 0;
        for (const auto& rd : ragged_dimensions_) {
            total_size += rd;
        }

        for (size_t i = N - 1; i-- > 0;) {
            strides_[i] = total_size;
            total_size *= dimensions_[i];
        }

        ragged_offsets_.resize(ragged_dimensions_.size() + 1, 0);
        for (size_t i = 0; i < ragged_dimensions_.size(); ++i) {
            ragged_offsets_.at(i + 1) = ragged_offsets_.at(i) + ragged_dimensions_.at(i);
        }

        data_.resize(total_size, default_value);
    }

    void inner_fill(const std::array<size_t, N - 1>& indices, const T& value) {
        const auto [begin, end] = inner_iterators(indices);
        std::fill(begin, end, value);
    }

    void inner_fill(const std::array<size_t, N - 1>& indices, const std::span<T const> values) {
#ifndef NDEBUG
        if (values.size() != ragged_dimensions_.at(indices.back())) {
            throw std::invalid_argument("Incorrect number of values");
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

private:
    std::vector<T> data_;
    std::vector<size_t> ragged_dimensions_;
    std::vector<size_t> ragged_offsets_;
    std::array<size_t, N - 1> dimensions_;
    std::array<size_t, N - 1> strides_;

    /// Helper function to calculate the 1D index from multi-dimensional indices
    /// @param indices The multi-dimensional indices.
    /// @return The calculated 1D index.
    /// @throws std::out_of_range if any index is out of bounds.
    inline size_t calculate_index(const std::array<size_t, N>& indices) const {
#ifndef NDEBUG
        for (size_t i = 0; i < indices.size() - 1; ++i) {
            if (indices[i] >= dimensions_[i]) {
                throw std::out_of_range("Index out of bounds");
            }
        }
        if (indices.back() >= ragged_dimensions_[indices[indices.size() - 1]]) {
            throw std::out_of_range("Index out of bounds");
        }
#endif
        size_t index = 0;
        for (size_t i = 0; i < indices.size() - 1; ++i) {
            index += indices[i] * strides_[i];
        }
        index += ragged_offsets_[indices[indices.size() - 1]] + indices.back();
        return index;
    }

    /// Helper function to calculate start and end indices for iterators
    /// @param outer_indices The indices of the outer dimensions.
    /// @return A pair of start and end indices for the iterators.
    std::pair<size_t, size_t> calculate_start_end_indices(const std::array<size_t, N - 1>& outer_indices) const {
        std::array<size_t, N> full_indices{0};
        std::copy(outer_indices.begin(), outer_indices.end(), full_indices.begin());
        const size_t start_index = calculate_index(full_indices);
        const size_t end_index = start_index + ragged_dimensions_.at(outer_indices.back());
        // std::cout << "Start index: " << start_index << std::endl;
        // std::cout << "End index: " << end_index << std::endl;
        // std::cout << "Ragged dimensions: ";
        // for (const auto &e : ragged_dimensions_)
        // {
        //     std::cout << e << " ";
        // }
        // std::cout << std::endl;
        // std::cout << "Ragged dimension: " << ragged_dimensions_.at(outer_indices.back()) << std::endl;
        // std::cout << "Outer indices: ";
        // for (const auto &e : outer_indices)
        // {
        //     std::cout << e << " ";
        // }
        // std::cout << std::endl;
        return {start_index, end_index};
    }
};
