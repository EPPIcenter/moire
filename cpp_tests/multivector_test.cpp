#include <iostream>
#include "multivector.h"
#include <vector>
#include <cassert>

void test_multivector_construction() {
    std::array<size_t, 3> dimensions = {2, 3, 4};
    MultiVector<int, 3> mv(dimensions);
    assert(mv.total_size() == 24);
    std::cout << "Construction test passed.\n";
}

void test_multivector_element_access() {
    std::array<size_t, 3> dimensions = {2, 3, 4};
    MultiVector<int, 3> mv(dimensions);
    mv.at({0, 0, 0}) = 42;
    mv.at({0, 0, 1}) = 43;
    mv.at({0, 0, 2}) = 44;
    mv.at({0, 0, 3}) = 45;
    mv.at({1, 2, 3}) = 46;
    assert(mv.at({0, 0, 0}) == 42);
    assert(mv.at({0, 0, 1}) == 43);
    assert(mv.at({0, 0, 2}) == 44);
    assert(mv.at({0, 0, 3}) == 45);
    assert(mv.at({1, 2, 3}) == 46);
    std::cout << "Element access test passed.\n";
}


void test_multivector_dimension_size() {
    std::array<size_t, 3> dimensions = {2, 3, 4};
    MultiVector<int, 3> mv(dimensions);
    assert(mv.size(0) == 2);
    assert(mv.size(1) == 3);
    assert(mv.size(2) == 4);
    std::cout << "Dimension size test passed.\n";
}

void test_multivector_iteration() {
    std::array<size_t, 3> dimensions = {2, 3, 4};
    MultiVector<int, 3> mv(dimensions);
    for (std::size_t i = 0; i < 4; ++i) {
        mv.at({0, 0, i}) = i;
    }
    int index = 0;
    for (auto it = mv.inner_begin({0, 0}); it != mv.inner_end({0, 0}); ++it) {
        assert(*it == index++);
    }
    std::cout << "Iteration test passed.\n";
}

// Test case for RaggedMultiVector construction
void test_ragged_multivector_construction() {
    std::array<size_t, 2> dimensions = {2, 3};
    std::vector<size_t> ragged_dimensions = {2, 3, 4};
    RaggedMultiVector<int, 3> raggedMultiVector(dimensions, ragged_dimensions);
    RaggedMultiVector<int, 1> raggedMultiVector2({}, ragged_dimensions);
    
    // Check if the total size is correct
    size_t expected_total_size = 2 * 3 * (2 + 3 + 4);
    assert(raggedMultiVector.total_size() == expected_total_size);
    assert(raggedMultiVector2.total_size() == 9);

    assert(raggedMultiVector.dimensions() == (std::array<size_t, 3>{2, 3, 3}));
    assert(raggedMultiVector2.dimensions() == (std::array<size_t, 1>{3}));

    std::cout << "RaggedMultiVector construction test passed." << std::endl;
}

// Test case for RaggedMultiVector element access
void test_ragged_multivector_element_access() {
    std::array<size_t, 2> dimensions = {2, 3};
    std::vector<size_t> ragged_dimensions = {2, 3, 4};
    RaggedMultiVector<int, 3> raggedMultiVector(dimensions, ragged_dimensions);
    RaggedMultiVector<float, 1> raggedMultiVector2({}, ragged_dimensions);

    raggedMultiVector.at({0, 0, 0, 0}) = 10;
    raggedMultiVector.at({0, 1, 1, 2}) = 20;
    raggedMultiVector.at({1, 2, 2, 3}) = 30;
    raggedMultiVector.at({0, 2, 0, 1}) = 40;
    raggedMultiVector.at({1, 0, 1, 0}) = 50;
    raggedMultiVector.at({1, 1, 2, 2}) = 60;

    raggedMultiVector2.at({0, 1}) = 10.0f;
    raggedMultiVector2.at({1, 2}) = 20.0f;
    raggedMultiVector2.at({2, 3}) = 30.0f; 

    assert(raggedMultiVector.at({0, 0, 0, 0}) == 10);
    assert(raggedMultiVector.at({0, 1, 1, 2}) == 20);
    assert(raggedMultiVector.at({1, 2, 2, 3}) == 30);
    assert(raggedMultiVector.at({0, 2, 0, 1}) == 40);
    assert(raggedMultiVector.at({1, 0, 1, 0}) == 50);
    assert(raggedMultiVector.at({1, 1, 2, 2}) == 60);

    assert(raggedMultiVector2.at({0, 1}) == 10.0f);
    assert(raggedMultiVector2.at({1, 2}) == 20.0f);
    assert(raggedMultiVector2.at({2, 3}) == 30.0f);

    assert(raggedMultiVector.total_size() == 54);
    assert(raggedMultiVector2.total_size() == 9);

    std::cout << "RaggedMultiVector element access test passed." << std::endl;
}

// Test case for RaggedMultiVector iteration
void test_ragged_multivector_iteration() {
    std::array<size_t, 2> dimensions = {2, 3};
    std::vector<size_t> ragged_dimensions = {2, 3, 4};
    RaggedMultiVector<int, 3> raggedMultiVector(dimensions, ragged_dimensions);

    for (size_t i = 0; i < ragged_dimensions[0]; ++i) {
        raggedMultiVector.at({0, 0, 0, i}) = i;
    }

    int index = 0;
    for (auto it = raggedMultiVector.inner_begin({0, 0, 0}); it != raggedMultiVector.inner_end({0, 0, 0}); ++it) {
        assert(*it == index++);
    }

    std::cout << "RaggedMultiVector iteration test passed." << std::endl;
}

int main() {
    test_multivector_construction();
    test_multivector_element_access();
    test_multivector_dimension_size();
    test_multivector_iteration();
    test_ragged_multivector_construction();
    test_ragged_multivector_element_access();
    test_ragged_multivector_iteration();
    std::cout << "All tests passed." << std::endl;
    return 0;
}


