#include "test_framework.hpp"
#include "multivector.h"
#include <vector>
#include <array>
#include <random>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <span>

namespace moire_test {

// Test suite for MultiVector
class MultiVectorTestSuite : public TestSuite {
public:
    MultiVectorTestSuite() : TestSuite("MultiVector") {
        add_test("Construction", [this]() { test_construction(); });
        add_test("Element Access", [this]() { test_element_access(); });
        add_test("Dimension Size", [this]() { test_dimensions(); });
        add_test("Inner Fill", [this]() { test_inner_fill(); });
        add_test("Data Access", [this]() { test_data_access(); });
        add_test("Copy Constructor", [this]() { test_copy_constructor(); });
        add_test("Assignment Operator", [this]() { test_assignment_operator(); });
        add_test("Move Constructor", [this]() { test_move_constructor(); });
        add_test("Move Assignment", [this]() { test_move_assignment(); });
        add_test("Large Dimensions", [this]() { test_large_dimensions(); });
        add_test("Edge Cases", [this]() { test_edge_cases(); });
        add_test("Performance", [this]() { test_performance(); });
        
        // Additional tests from old framework
        add_test("Iteration", [this]() { test_iteration(); });
        add_test("Reduction Operations", [this]() { test_reduction_operations(); });
        add_test("Logsumexp Reduction", [this]() { test_logsumexp_reduction(); });
        add_test("Transform Operations", [this]() { test_transform_operations(); });
        add_test("Element Transform Operations", [this]() { test_element_transform_operations(); });
        add_test("RaggedMultiVector Construction", [this]() { test_ragged_multivector_construction(); });
        add_test("RaggedMultiVector Element Access", [this]() { test_ragged_multivector_element_access(); });
        add_test("RaggedMultiVector Iteration", [this]() { test_ragged_multivector_iteration(); });
        add_test("N=1 Specialization", [this]() { test_n1_specialization(); });
        add_test("Elementwise Operations", [this]() { test_elementwise_operations(); });
        add_test("Softmax Operations", [this]() { test_softmax_operations(); });
        add_test("Invariant Default Constructor", [this]() { test_invariant_default_constructor(); });
        add_test("Invariant After Resize", [this]() { test_invariant_after_resize(); });
        add_test("Indexing N1 N2 N3", [this]() { test_indexing_n1_n2_n3(); });
        add_test("Parallel Sequential Equivalence", [this]() { test_parallel_sequential_equivalence(); });
        add_test("Logsumexp Numerical Stability", [this]() { test_logsumexp_numerical_stability(); });
        add_test("Softmax Normalization", [this]() { test_softmax_normalization(); });
        add_test("Ragged Bounds And Offsets", [this]() { test_ragged_bounds_and_offsets(); });
        add_test("At Bounds Check", [this]() { test_at_bounds_check(); });
    }

private:
    void test_construction() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv(dimensions);
        
        ASSERT_EQ(mv.total_size(), 24);
        auto dims = mv.dimensions();
        ASSERT_EQ(dims[0], 2);
        ASSERT_EQ(dims[1], 3);
        ASSERT_EQ(dims[2], 4);
    }

    void test_element_access() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv(dimensions);
        
        // Test setting and getting elements
        mv.at({0, 0, 0}) = 42;
        mv.at({0, 0, 1}) = 43;
        mv.at({0, 0, 2}) = 44;
        mv.at({0, 0, 3}) = 45;
        mv.at({1, 2, 3}) = 46;
        
        ASSERT_EQ(mv.at({0, 0, 0}), 42);
        ASSERT_EQ(mv.at({0, 0, 1}), 43);
        ASSERT_EQ(mv.at({0, 0, 2}), 44);
        ASSERT_EQ(mv.at({0, 0, 3}), 45);
        ASSERT_EQ(mv.at({1, 2, 3}), 46);
    }

    void test_dimensions() {
        std::array<size_t, 4> dimensions = {2, 3, 4, 5};
        MultiVector<int, 4> mv(dimensions);
        
        auto dims = mv.dimensions();
        ASSERT_EQ(dims[0], 2);
        ASSERT_EQ(dims[1], 3);
        ASSERT_EQ(dims[2], 4);
        ASSERT_EQ(dims[3], 5);
        ASSERT_EQ(mv.total_size(), 120);
    }

    void test_inner_fill() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv(dimensions);
        
        // Fill inner dimensions
        mv.inner_fill({0, 0}, 1);
        mv.inner_fill({0, 1}, 2);
        mv.inner_fill({0, 2}, 3);
        mv.inner_fill({1, 0}, 4);
        mv.inner_fill({1, 1}, 5);
        mv.inner_fill({1, 2}, 6);
        
        // Verify the fills
        ASSERT_EQ(mv.at({0, 0, 0}), 1);
        ASSERT_EQ(mv.at({0, 0, 1}), 1);
        ASSERT_EQ(mv.at({0, 0, 2}), 1);
        ASSERT_EQ(mv.at({0, 0, 3}), 1);
        ASSERT_EQ(mv.at({0, 1, 0}), 2);
        ASSERT_EQ(mv.at({0, 1, 1}), 2);
        ASSERT_EQ(mv.at({0, 1, 2}), 2);
        ASSERT_EQ(mv.at({0, 1, 3}), 2);
    }

    void test_data_access() {
        std::array<size_t, 2> dimensions = {3, 4};
        MultiVector<int, 2> mv(dimensions);
        
        // Fill with sequential values using at() method
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                mv.at({i, j}) = static_cast<int>(i * 4 + j);
            }
        }
        
        // Verify access through at()
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 4; ++j) {
                size_t expected = i * 4 + j;
                ASSERT_EQ(mv.at({i, j}), static_cast<int>(expected));
            }
        }
    }

    void test_copy_constructor() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv1(dimensions);
        
        // Fill with data using at() method
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    mv1.at({i, j, k}) = static_cast<int>(i * 12 + j * 4 + k);
                }
            }
        }
        
        // Copy construct
        MultiVector<int, 3> mv2(mv1);
        
        // Verify they're equal
        ASSERT_EQ(mv1.total_size(), mv2.total_size());
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    ASSERT_EQ(mv1.at({i, j, k}), mv2.at({i, j, k}));
                }
            }
        }
        
        // Verify they're independent
        mv1.at({0, 0, 0}) = 999;
        ASSERT_NE(mv1.at({0, 0, 0}), mv2.at({0, 0, 0}));
    }

    void test_assignment_operator() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv1(dimensions);
        MultiVector<int, 3> mv2(dimensions);
        
        // Fill mv1 with data using at() method
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    mv1.at({i, j, k}) = static_cast<int>(i * 12 + j * 4 + k);
                }
            }
        }
        
        // Assign mv1 to mv2
        mv2 = mv1;
        
        // Verify they're equal
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    ASSERT_EQ(mv1.at({i, j, k}), mv2.at({i, j, k}));
                }
            }
        }
        
        // Verify they're independent
        mv1.at({0, 0, 0}) = 999;
        ASSERT_NE(mv1.at({0, 0, 0}), mv2.at({0, 0, 0}));
    }

    void test_move_constructor() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv1(dimensions);
        
        // Fill with data using at() method
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    mv1.at({i, j, k}) = static_cast<int>(i * 12 + j * 4 + k);
                }
            }
        }
        
        // Move construct
        MultiVector<int, 3> mv2(std::move(mv1));
        
        // Verify mv2 has the data
        ASSERT_EQ(mv2.total_size(), 24);
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    ASSERT_EQ(mv2.at({i, j, k}), static_cast<int>(i * 12 + j * 4 + k));
                }
            }
        }
    }

    void test_move_assignment() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv1(dimensions);
        MultiVector<int, 3> mv2(dimensions);
        
        // Fill mv1 with data using at() method
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    mv1.at({i, j, k}) = static_cast<int>(i * 12 + j * 4 + k);
                }
            }
        }
        
        // Move assign
        mv2 = std::move(mv1);
        
        // Verify mv2 has the data
        for (size_t i = 0; i < 2; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < 4; ++k) {
                    ASSERT_EQ(mv2.at({i, j, k}), static_cast<int>(i * 12 + j * 4 + k));
                }
            }
        }
    }

    void test_large_dimensions() {
        std::array<size_t, 3> dimensions = {100, 100, 100};
        MultiVector<int, 3> mv(dimensions);
        
        ASSERT_EQ(mv.total_size(), 1000000);
        auto dims = mv.dimensions();
        ASSERT_EQ(dims[0], 100);
        ASSERT_EQ(dims[1], 100);
        ASSERT_EQ(dims[2], 100);
        
        // Test access to some elements
        mv.at({0, 0, 0}) = 42;
        mv.at({99, 99, 99}) = 43;
        ASSERT_EQ(mv.at({0, 0, 0}), 42);
        ASSERT_EQ(mv.at({99, 99, 99}), 43);
    }

    void test_edge_cases() {
        // Test with dimension size 1
        std::array<size_t, 3> dimensions = {1, 1, 1};
        MultiVector<int, 3> mv1(dimensions);
        ASSERT_EQ(mv1.total_size(), 1);
        mv1.at({0, 0, 0}) = 42;
        ASSERT_EQ(mv1.at({0, 0, 0}), 42);
        
        // Test with very large single dimension
        std::array<size_t, 1> large_dim = {1000000};
        MultiVector<int, 1> mv2(large_dim);
        ASSERT_EQ(mv2.total_size(), 1000000);
        mv2.at({0}) = 1;
        mv2.at({999999}) = 2;
        ASSERT_EQ(mv2.at({0}), 1);
        ASSERT_EQ(mv2.at({999999}), 2);
    }

    void test_performance() {
        std::array<size_t, 3> dimensions = {100, 100, 100};
        MultiVector<int, 3> mv(dimensions);
        
        // Test write performance
        auto start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < 1000; ++i) {
            for (size_t j = 0; j < 1000; ++j) {
                for (size_t k = 0; k < 1000; ++k) {
                    if (i < 100 && j < 100 && k < 100) {
                        mv.at({i, j, k}) = static_cast<int>(i + j + k);
                    }
                }
            }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        // Performance should be reasonable (less than 1 second for this test)
        ASSERT_LT(duration, 1000);
    }

    void test_iteration() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv(dimensions);
        
        // Fill inner dimension with sequential values
        for (std::size_t i = 0; i < 4; ++i) {
            mv.at({0, 0, i}) = static_cast<int>(i);
        }
        
        // Test iteration over inner dimension
        int index = 0;
        for (auto it = mv.inner_begin({0, 0}); it != mv.inner_end({0, 0}); ++it) {
            ASSERT_EQ(*it, index++);
        }
    }

    void test_reduction_operations() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv(dimensions);
        
        // Fill with test data
        // First slice (i=0)
        mv.at({0, 0, 0}) = 1;  mv.at({0, 0, 1}) = 2;  mv.at({0, 0, 2}) = 3;  mv.at({0, 0, 3}) = 4;
        mv.at({0, 1, 0}) = 2;  mv.at({0, 1, 1}) = 3;  mv.at({0, 1, 2}) = 4;  mv.at({0, 1, 3}) = 5;
        mv.at({0, 2, 0}) = 3;  mv.at({0, 2, 1}) = 4;  mv.at({0, 2, 2}) = 5;  mv.at({0, 2, 3}) = 6;
        
        // Second slice (i=1)
        mv.at({1, 0, 0}) = 4;  mv.at({1, 0, 1}) = 5;  mv.at({1, 0, 2}) = 6;  mv.at({1, 0, 3}) = 7;
        mv.at({1, 1, 0}) = 5;  mv.at({1, 1, 1}) = 6;  mv.at({1, 1, 2}) = 7;  mv.at({1, 1, 3}) = 8;
        mv.at({1, 2, 0}) = 6;  mv.at({1, 2, 1}) = 7;  mv.at({1, 2, 2}) = 8;  mv.at({1, 2, 3}) = 9;

        // Test sum reduction
        auto sums = mv.sum();
        ASSERT_EQ(sums.dimensions(), (std::array<size_t, 2>{2, 3}));
        ASSERT_EQ(sums.at({0, 0}), 10);  // 1+2+3+4
        ASSERT_EQ(sums.at({0, 1}), 14);  // 2+3+4+5
        ASSERT_EQ(sums.at({0, 2}), 18);  // 3+4+5+6
        ASSERT_EQ(sums.at({1, 0}), 22);  // 4+5+6+7
        ASSERT_EQ(sums.at({1, 1}), 26);  // 5+6+7+8
        ASSERT_EQ(sums.at({1, 2}), 30);  // 6+7+8+9

        // Test product reduction
        auto products = mv.product();
        ASSERT_EQ(products.dimensions(), (std::array<size_t, 2>{2, 3}));
        ASSERT_EQ(products.at({0, 0}), 24);   // 1*2*3*4
        ASSERT_EQ(products.at({0, 1}), 120);  // 2*3*4*5
        ASSERT_EQ(products.at({0, 2}), 360);  // 3*4*5*6
        ASSERT_EQ(products.at({1, 0}), 840);  // 4*5*6*7
        ASSERT_EQ(products.at({1, 1}), 1680); // 5*6*7*8
        ASSERT_EQ(products.at({1, 2}), 3024); // 6*7*8*9

        // Test max reduction
        auto maxes = mv.max();
        ASSERT_EQ(maxes.dimensions(), (std::array<size_t, 2>{2, 3}));
        ASSERT_EQ(maxes.at({0, 0}), 4);  // max(1,2,3,4)
        ASSERT_EQ(maxes.at({0, 1}), 5);  // max(2,3,4,5)
        ASSERT_EQ(maxes.at({0, 2}), 6);  // max(3,4,5,6)
        ASSERT_EQ(maxes.at({1, 0}), 7);  // max(4,5,6,7)
        ASSERT_EQ(maxes.at({1, 1}), 8);  // max(5,6,7,8)
        ASSERT_EQ(maxes.at({1, 2}), 9);  // max(6,7,8,9)

        // Test min reduction
        auto mins = mv.min();
        ASSERT_EQ(mins.dimensions(), (std::array<size_t, 2>{2, 3}));
        ASSERT_EQ(mins.at({0, 0}), 1);  // min(1,2,3,4)
        ASSERT_EQ(mins.at({0, 1}), 2);  // min(2,3,4,5)
        ASSERT_EQ(mins.at({0, 2}), 3);  // min(3,4,5,6)
        ASSERT_EQ(mins.at({1, 0}), 4);  // min(4,5,6,7)
        ASSERT_EQ(mins.at({1, 1}), 5);  // min(5,6,7,8)
        ASSERT_EQ(mins.at({1, 2}), 6);  // min(6,7,8,9)
    }

    void test_logsumexp_reduction() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<double, 3> mv_double(dimensions);
        
        // Fill with test data
        // First slice (i=0)
        mv_double.at({0, 0, 0}) = 1.0;  mv_double.at({0, 0, 1}) = 2.0;  mv_double.at({0, 0, 2}) = 3.0;  mv_double.at({0, 0, 3}) = 4.0;
        mv_double.at({0, 1, 0}) = 2.0;  mv_double.at({0, 1, 1}) = 3.0;  mv_double.at({0, 1, 2}) = 4.0;  mv_double.at({0, 1, 3}) = 5.0;
        mv_double.at({0, 2, 0}) = 3.0;  mv_double.at({0, 2, 1}) = 4.0;  mv_double.at({0, 2, 2}) = 5.0;  mv_double.at({0, 2, 3}) = 6.0;
        
        // Second slice (i=1)
        mv_double.at({1, 0, 0}) = 4.0;  mv_double.at({1, 0, 1}) = 5.0;  mv_double.at({1, 0, 2}) = 6.0;  mv_double.at({1, 0, 3}) = 7.0;
        mv_double.at({1, 1, 0}) = 5.0;  mv_double.at({1, 1, 1}) = 6.0;  mv_double.at({1, 1, 2}) = 7.0;  mv_double.at({1, 1, 3}) = 8.0;
        mv_double.at({1, 2, 0}) = 6.0;  mv_double.at({1, 2, 1}) = 7.0;  mv_double.at({1, 2, 2}) = 8.0;  mv_double.at({1, 2, 3}) = 9.0;
        
        // Test logsumexp reduction
        auto logsumexps = mv_double.logsumexp();
        ASSERT_EQ(logsumexps.dimensions(), (std::array<size_t, 2>{2, 3}));
        
        // Calculate expected values manually
        // For {0, 0}: log(exp(1) + exp(2) + exp(3) + exp(4)) = log(e + e^2 + e^3 + e^4) ≈ 4.4402
        ASSERT_LT(std::abs(logsumexps.at({0, 0}) - 4.4402), 0.0001);
        
        // For {0, 1}: log(exp(2) + exp(3) + exp(4) + exp(5)) = log(e^2 + e^3 + e^4 + e^5) ≈ 5.4402
        ASSERT_LT(std::abs(logsumexps.at({0, 1}) - 5.4402), 0.0001);
        
        // Test parallel_logsumexp
        auto parallel_logsumexps = mv_double.parallel_logsumexp();
        ASSERT_EQ(parallel_logsumexps.dimensions(), (std::array<size_t, 2>{2, 3}));
        
        // Verify that parallel results match sequential results
        ASSERT_LT(std::abs(parallel_logsumexps.at({0, 0}) - logsumexps.at({0, 0})), 0.0001);
        ASSERT_LT(std::abs(parallel_logsumexps.at({0, 1}) - logsumexps.at({0, 1})), 0.0001);
    }

    void test_transform_operations() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv(dimensions);
        
        // Fill with test data
        mv.at({0, 0, 0}) = 1;  mv.at({0, 0, 1}) = 2;  mv.at({0, 0, 2}) = 3;  mv.at({0, 0, 3}) = 4;
        mv.at({0, 1, 0}) = 2;  mv.at({0, 1, 1}) = 3;  mv.at({0, 1, 2}) = 4;  mv.at({0, 1, 3}) = 5;
        mv.at({0, 2, 0}) = 3;  mv.at({0, 2, 1}) = 4;  mv.at({0, 2, 2}) = 5;  mv.at({0, 2, 3}) = 6;
        mv.at({1, 0, 0}) = 4;  mv.at({1, 0, 1}) = 5;  mv.at({1, 0, 2}) = 6;  mv.at({1, 0, 3}) = 7;
        mv.at({1, 1, 0}) = 5;  mv.at({1, 1, 1}) = 6;  mv.at({1, 1, 2}) = 7;  mv.at({1, 1, 3}) = 8;
        mv.at({1, 2, 0}) = 6;  mv.at({1, 2, 1}) = 7;  mv.at({1, 2, 2}) = 8;  mv.at({1, 2, 3}) = 9;
        
        // Create a copy for testing
        MultiVector<int, 3> mv_copy = mv;
        
        // Create values for each inner slice (2*3 = 6 values)
        std::vector<int> add_values = {1, 2, 3, 4, 5, 6};
        
        // Test add operation
        mv = mv.add(add_values);
        
        // Verify results
        ASSERT_EQ(mv.at({0, 0, 0}), 2);  // 1+1
        ASSERT_EQ(mv.at({0, 0, 1}), 3);  // 2+1
        ASSERT_EQ(mv.at({0, 0, 2}), 4);  // 3+1
        ASSERT_EQ(mv.at({0, 0, 3}), 5);  // 4+1
        
        // Test subtract operation
        mv = mv_copy;
        mv = mv.subtract(add_values);
        
        // Verify results
        ASSERT_EQ(mv.at({0, 0, 0}), 0);  // 1-1
        ASSERT_EQ(mv.at({0, 0, 1}), 1);  // 2-1
        ASSERT_EQ(mv.at({0, 0, 2}), 2);  // 3-1
        ASSERT_EQ(mv.at({0, 0, 3}), 3);  // 4-1
        
        // Test multiply operation
        mv = mv_copy;
        std::vector<int> mult_values = {2, 3, 4, 5, 6, 7};
        mv = mv.multiply(mult_values);
        
        // Verify results
        ASSERT_EQ(mv.at({0, 0, 0}), 2);   // 1*2
        ASSERT_EQ(mv.at({0, 0, 1}), 4);   // 2*2
        ASSERT_EQ(mv.at({0, 0, 2}), 6);   // 3*2
        ASSERT_EQ(mv.at({0, 0, 3}), 8);   // 4*2
        
        // Test divide operation
        mv = mv_copy;
        std::vector<int> div_values = {2, 2, 2, 2, 2, 2};
        mv = mv.divide(div_values);
        
        // Verify results (integer division)
        ASSERT_EQ(mv.at({0, 0, 0}), 0);  // 1/2
        ASSERT_EQ(mv.at({0, 0, 1}), 1);  // 2/2
        ASSERT_EQ(mv.at({0, 0, 2}), 1);  // 3/2
        ASSERT_EQ(mv.at({0, 0, 3}), 2);  // 4/2
    }

    void test_element_transform_operations() {
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv(dimensions);
        
        // Fill with test data
        mv.at({0, 0, 0}) = 1;  mv.at({0, 0, 1}) = 2;  mv.at({0, 0, 2}) = 3;  mv.at({0, 0, 3}) = 4;
        mv.at({0, 1, 0}) = 2;  mv.at({0, 1, 1}) = 3;  mv.at({0, 1, 2}) = 4;  mv.at({0, 1, 3}) = 5;
        mv.at({0, 2, 0}) = 3;  mv.at({0, 2, 1}) = 4;  mv.at({0, 2, 2}) = 5;  mv.at({0, 2, 3}) = 6;
        mv.at({1, 0, 0}) = 4;  mv.at({1, 0, 1}) = 5;  mv.at({1, 0, 2}) = 6;  mv.at({1, 0, 3}) = 7;
        mv.at({1, 1, 0}) = 5;  mv.at({1, 1, 1}) = 6;  mv.at({1, 1, 2}) = 7;  mv.at({1, 1, 3}) = 8;
        mv.at({1, 2, 0}) = 6;  mv.at({1, 2, 1}) = 7;  mv.at({1, 2, 2}) = 8;  mv.at({1, 2, 3}) = 9;
        
        // Create a copy for testing
        MultiVector<int, 3> mv_copy = mv;
        
        // Create values for each element position in the inner dimension (4 values)
        std::vector<int> add_values = {1, 2, 3, 4};
        
        // Test element_add operation
        mv = mv.element_add(add_values);
        
        // Verify results
        ASSERT_EQ(mv.at({0, 0, 0}), 2);  // 1+1
        ASSERT_EQ(mv.at({0, 0, 1}), 4);  // 2+2
        ASSERT_EQ(mv.at({0, 0, 2}), 6);  // 3+3
        ASSERT_EQ(mv.at({0, 0, 3}), 8);  // 4+4
        
        // Test element_subtract operation
        mv = mv_copy;
        mv = mv.element_subtract(add_values);
        
        // Verify results
        ASSERT_EQ(mv.at({0, 0, 0}), 0);  // 1-1
        ASSERT_EQ(mv.at({0, 0, 1}), 0);  // 2-2
        ASSERT_EQ(mv.at({0, 0, 2}), 0);  // 3-3
        ASSERT_EQ(mv.at({0, 0, 3}), 0);  // 4-4
        
        // Test element_multiply operation
        mv = mv_copy;
        std::vector<int> mult_values = {2, 3, 4, 5};
        mv = mv.element_multiply(mult_values);
        
        // Verify results
        ASSERT_EQ(mv.at({0, 0, 0}), 2);   // 1*2
        ASSERT_EQ(mv.at({0, 0, 1}), 6);   // 2*3
        ASSERT_EQ(mv.at({0, 0, 2}), 12);  // 3*4
        ASSERT_EQ(mv.at({0, 0, 3}), 20);  // 4*5
    }

    void test_ragged_multivector_construction() {
        std::array<size_t, 3> dimensions = {2, 3, 3};
        std::vector<size_t> ragged_dimensions = {2, 3, 4};
        RaggedMultiVector<int, 4> raggedMultiVector(dimensions, std::span<size_t const>(ragged_dimensions));
        
        // Check if the total size is correct
        size_t expected_total_size = 2 * 3 * (2 + 3 + 4);
        ASSERT_EQ(raggedMultiVector.total_size(), expected_total_size);

        // Check dimensions
        auto dims = raggedMultiVector.dimensions();
        ASSERT_EQ(dims[0], 2);
        ASSERT_EQ(dims[1], 3);
        ASSERT_EQ(dims[2], 3); // The last dimension is the size of ragged_dimensions
    }

    void test_ragged_multivector_element_access() {
        std::array<size_t, 2> dimensions = {2, 3}; // 2 samples, 3 loci, for example
        std::vector<size_t> ragged_dimensions = {2, 3, 4};
        RaggedMultiVector<int, 3> raggedMultiVector(dimensions, std::span<size_t const>(ragged_dimensions));

        // Access elements using the at method
        raggedMultiVector.at({0, 0, 0}) = 10;
        raggedMultiVector.at({0, 1, 1}) = 20;
        raggedMultiVector.at({1, 2, 2}) = 30;
        raggedMultiVector.at({0, 2, 0}) = 40;
        raggedMultiVector.at({1, 0, 1}) = 50;
        raggedMultiVector.at({1, 1, 2}) = 60;

        // Verify the values
        ASSERT_EQ(raggedMultiVector.at({0, 0, 0}), 10);
        ASSERT_EQ(raggedMultiVector.at({0, 1, 1}), 20);
        ASSERT_EQ(raggedMultiVector.at({1, 2, 2}), 30);
        ASSERT_EQ(raggedMultiVector.at({0, 2, 0}), 40);
        ASSERT_EQ(raggedMultiVector.at({1, 0, 1}), 50);
        ASSERT_EQ(raggedMultiVector.at({1, 1, 2}), 60);

        // Check total size
        ASSERT_EQ(raggedMultiVector.total_size(), 18);
    }

    void test_ragged_multivector_iteration() {
        std::array<size_t, 2> dimensions = {2, 3};
        std::vector<size_t> ragged_dimensions = {2, 3, 4};
        RaggedMultiVector<int, 3> raggedMultiVector(dimensions, std::span<size_t const>(ragged_dimensions));

        // Fill the first inner dimension with values
        for (size_t i = 0; i < ragged_dimensions[0]; ++i) {
            raggedMultiVector.at({0, 0, i}) = static_cast<int>(i);
        }

        // Iterate over the inner dimension and verify values
        int index = 0;
        for (auto it = raggedMultiVector.inner_begin({0, 0}); it != raggedMultiVector.inner_end({0, 0}); ++it) {
            ASSERT_EQ(*it, index++);
        }
    }

    void test_n1_specialization() {
        // Create a 1D MultiVector with 5 elements
        std::array<size_t, 1> dimensions = {5};
        MultiVector<int, 1> mv(dimensions);
        
        // Fill with test data
        mv.at({0}) = 1;
        mv.at({1}) = 2;
        mv.at({2}) = 3;
        mv.at({3}) = 4;
        mv.at({4}) = 5;
        
        // Test sum reduction - should return a scalar
        int sum = mv.sum();
        ASSERT_EQ(sum, 15);  // 1+2+3+4+5
        
        // Test product reduction - should return a scalar
        int product = mv.product();
        ASSERT_EQ(product, 120);  // 1*2*3*4*5
        
        // Test max reduction - should return a scalar
        int max_val = mv.max();
        ASSERT_EQ(max_val, 5);  // max(1,2,3,4,5)
        
        // Test min reduction - should return a scalar
        int min_val = mv.min();
        ASSERT_EQ(min_val, 1);  // min(1,2,3,4,5)
        
        // Test parallel reductions
        int parallel_sum = mv.parallel_sum();
        ASSERT_EQ(parallel_sum, 15);
        
        int parallel_product = mv.parallel_product();
        ASSERT_EQ(parallel_product, 120);
        
        int parallel_max = mv.parallel_max();
        ASSERT_EQ(parallel_max, 5);
        
        int parallel_min = mv.parallel_min();
        ASSERT_EQ(parallel_min, 1);
    }

    void test_elementwise_operations() {
        // Create two 2x3x4 MultiVectors
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<int, 3> mv1(dimensions);
        MultiVector<int, 3> mv2(dimensions);
        
        // Fill with test data
        mv1.at({0, 0, 0}) = 1;  mv1.at({0, 0, 1}) = 2;  mv1.at({0, 0, 2}) = 3;  mv1.at({0, 0, 3}) = 4;
        mv1.at({0, 1, 0}) = 2;  mv1.at({0, 1, 1}) = 3;  mv1.at({0, 1, 2}) = 4;  mv1.at({0, 1, 3}) = 5;
        mv1.at({0, 2, 0}) = 3;  mv1.at({0, 2, 1}) = 4;  mv1.at({0, 2, 2}) = 5;  mv1.at({0, 2, 3}) = 6;
        mv1.at({1, 0, 0}) = 4;  mv1.at({1, 0, 1}) = 5;  mv1.at({1, 0, 2}) = 6;  mv1.at({1, 0, 3}) = 7;
        mv1.at({1, 1, 0}) = 5;  mv1.at({1, 1, 1}) = 6;  mv1.at({1, 1, 2}) = 7;  mv1.at({1, 1, 3}) = 8;
        mv1.at({1, 2, 0}) = 6;  mv1.at({1, 2, 1}) = 7;  mv1.at({1, 2, 2}) = 8;  mv1.at({1, 2, 3}) = 9;
        
        // Fill second multivector with different values
        mv2.at({0, 0, 0}) = 10;  mv2.at({0, 0, 1}) = 11;  mv2.at({0, 0, 2}) = 12;  mv2.at({0, 0, 3}) = 13;
        mv2.at({0, 1, 0}) = 11;  mv2.at({0, 1, 1}) = 12;  mv2.at({0, 1, 2}) = 13;  mv2.at({0, 1, 3}) = 14;
        mv2.at({0, 2, 0}) = 12;  mv2.at({0, 2, 1}) = 13;  mv2.at({0, 2, 2}) = 14;  mv2.at({0, 2, 3}) = 15;
        mv2.at({1, 0, 0}) = 13;  mv2.at({1, 0, 1}) = 14;  mv2.at({1, 0, 2}) = 15;  mv2.at({1, 0, 3}) = 16;
        mv2.at({1, 1, 0}) = 14;  mv2.at({1, 1, 1}) = 15;  mv2.at({1, 1, 2}) = 16;  mv2.at({1, 1, 3}) = 17;
        mv2.at({1, 2, 0}) = 15;  mv2.at({1, 2, 1}) = 16;  mv2.at({1, 2, 2}) = 17;  mv2.at({1, 2, 3}) = 18;
        
        // Test addition
        MultiVector<int, 3> sum = mv1 + mv2;
        
        // Verify results
        ASSERT_EQ(sum.at({0, 0, 0}), 11);  // 1+10
        ASSERT_EQ(sum.at({0, 0, 1}), 13);  // 2+11
        ASSERT_EQ(sum.at({0, 0, 2}), 15);  // 3+12
        ASSERT_EQ(sum.at({0, 0, 3}), 17);  // 4+13
        
        // Test subtraction
        MultiVector<int, 3> diff = mv2 - mv1;
        
        // Verify results
        ASSERT_EQ(diff.at({0, 0, 0}), 9);   // 10-1
        ASSERT_EQ(diff.at({0, 0, 1}), 9);   // 11-2
        ASSERT_EQ(diff.at({0, 0, 2}), 9);   // 12-3
        ASSERT_EQ(diff.at({0, 0, 3}), 9);   // 13-4
        
        // Test multiplication
        MultiVector<int, 3> product = mv1 * mv2;
        
        // Verify results
        ASSERT_EQ(product.at({0, 0, 0}), 10);   // 1*10
        ASSERT_EQ(product.at({0, 0, 1}), 22);   // 2*11
        ASSERT_EQ(product.at({0, 0, 2}), 36);   // 3*12
        ASSERT_EQ(product.at({0, 0, 3}), 52);   // 4*13
        
        // Test division (integer division)
        MultiVector<int, 3> quotient = mv2 / mv1;
        
        // Verify results
        ASSERT_EQ(quotient.at({0, 0, 0}), 10);  // 10/1
        ASSERT_EQ(quotient.at({0, 0, 1}), 5);   // 11/2
        ASSERT_EQ(quotient.at({0, 0, 2}), 4);   // 12/3
        ASSERT_EQ(quotient.at({0, 0, 3}), 3);   // 13/4
    }

    void test_softmax_operations() {
        // Create a 2x3x4 MultiVector with floating point values
        std::array<size_t, 3> dimensions = {2, 3, 4};
        MultiVector<double, 3> mv(dimensions);
        
        // Fill with test data
        mv.at({0, 0, 0}) = 1.0;  mv.at({0, 0, 1}) = 2.0;  mv.at({0, 0, 2}) = 3.0;  mv.at({0, 0, 3}) = 4.0;
        mv.at({0, 1, 0}) = 2.0;  mv.at({0, 1, 1}) = 3.0;  mv.at({0, 1, 2}) = 4.0;  mv.at({0, 1, 3}) = 5.0;
        mv.at({0, 2, 0}) = 3.0;  mv.at({0, 2, 1}) = 4.0;  mv.at({0, 2, 2}) = 5.0;  mv.at({0, 2, 3}) = 6.0;
        mv.at({1, 0, 0}) = 4.0;  mv.at({1, 0, 1}) = 5.0;  mv.at({1, 0, 2}) = 6.0;  mv.at({1, 0, 3}) = 7.0;
        mv.at({1, 1, 0}) = 5.0;  mv.at({1, 1, 1}) = 6.0;  mv.at({1, 1, 2}) = 7.0;  mv.at({1, 1, 3}) = 8.0;
        mv.at({1, 2, 0}) = 6.0;  mv.at({1, 2, 1}) = 7.0;  mv.at({1, 2, 2}) = 8.0;  mv.at({1, 2, 3}) = 9.0;
        
        // Test regular softmax
        auto softmax_result = mv.softmax(false);
        
        // Verify results for first inner slice
        double sum0 = 0.0;
        for (size_t k = 0; k < 4; ++k) {
            sum0 += softmax_result.at({0, 0, k});
        }
        ASSERT_LT(std::abs(sum0 - 1.0), 1e-10);  // Sum should be 1.0
        
        // Test log-space softmax
        auto log_softmax_result = mv.softmax(true);
        
        // Verify results for first inner slice
        double log_sum0 = 0.0;
        for (size_t k = 0; k < 4; ++k) {
            log_sum0 += std::exp(log_softmax_result.at({0, 0, k}));
        }
        ASSERT_LT(std::abs(log_sum0 - 1.0), 1e-10);  // Sum should be 1.0
        
        // Test parallel softmax
        auto parallel_softmax_result = mv.parallel_softmax(false);
        
        // Verify results for first inner slice
        double parallel_sum0 = 0.0;
        for (size_t k = 0; k < 4; ++k) {
            parallel_sum0 += parallel_softmax_result.at({0, 0, k});
        }
        ASSERT_LT(std::abs(parallel_sum0 - 1.0), 1e-10);  // Sum should be 1.0
    }

    void test_invariant_default_constructor() {
        MultiVector<int, 3> mv3;
        ASSERT_EQ(mv3.total_size(), 0u);
        auto d3 = mv3.dimensions();
        size_t product = 1;
        for (size_t i = 0; i < 3; ++i) product *= d3[i];
        ASSERT_EQ(product, 0u);

        MultiVector<double, 1> mv1;
        ASSERT_EQ(mv1.total_size(), 0u);
        ASSERT_EQ(mv1.dimensions()[0], 0u);

        RaggedMultiVector<int, 3> rmv;
        ASSERT_EQ(rmv.total_size(), 0u);
    }

    void test_invariant_after_resize() {
        std::array<size_t, 3> dims = {2, 3, 4};
        MultiVector<int, 3> mv(dims);
        ASSERT_EQ(mv.total_size(), 24u);
        mv.resize({5, 2, 3});
        ASSERT_EQ(mv.total_size(), 30u);
        auto d = mv.dimensions();
        ASSERT_EQ(d[0], 5u);
        ASSERT_EQ(d[1], 2u);
        ASSERT_EQ(d[2], 3u);
    }

    void test_indexing_n1_n2_n3() {
        std::array<size_t, 1> d1 = {7};
        MultiVector<int, 1> mv1(d1);
        for (size_t i = 0; i < 7; ++i) mv1.at({i}) = static_cast<int>(i);
        for (size_t i = 0; i < 7; ++i) ASSERT_EQ(mv1.at({i}), static_cast<int>(i));

        std::array<size_t, 2> d2 = {3, 4};
        MultiVector<int, 2> mv2(d2);
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 4; ++j)
                mv2.at({i, j}) = static_cast<int>(i * 4 + j);
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 4; ++j)
                ASSERT_EQ(mv2.at({i, j}), static_cast<int>(i * 4 + j));

        std::array<size_t, 3> d3 = {2, 3, 4};
        MultiVector<int, 3> mv3(d3);
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 3; ++j)
                for (size_t k = 0; k < 4; ++k)
                    mv3.at({i, j, k}) = static_cast<int>(i * 12 + j * 4 + k);
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 3; ++j)
                for (size_t k = 0; k < 4; ++k)
                    ASSERT_EQ(mv3.at({i, j, k}), static_cast<int>(i * 12 + j * 4 + k));
    }

    void test_parallel_sequential_equivalence() {
        const unsigned seed = 42;
        std::mt19937 gen(seed);
        std::uniform_real_distribution<double> dis(0.1, 2.0);
        std::array<size_t, 3> dims = {4, 5, 6};
        MultiVector<double, 3> mv(dims);
        for (size_t i = 0; i < dims[0]; ++i)
            for (size_t j = 0; j < dims[1]; ++j)
                for (size_t k = 0; k < dims[2]; ++k)
                    mv.at({i, j, k}) = dis(gen);

        auto seq_sum = mv.sum(std::execution::seq);
        auto par_sum = mv.parallel_sum();
        for (size_t i = 0; i < dims[0]; ++i)
            for (size_t j = 0; j < dims[1]; ++j)
                ASSERT_LT(std::abs(seq_sum.at({i, j}) - par_sum.at({i, j})), 1e-9);

        auto seq_logsumexp = mv.logsumexp();
        auto par_logsumexp = mv.parallel_logsumexp();
        for (size_t i = 0; i < dims[0]; ++i)
            for (size_t j = 0; j < dims[1]; ++j)
                ASSERT_LT(std::abs(seq_logsumexp.at({i, j}) - par_logsumexp.at({i, j})), 1e-9);
    }

    void test_logsumexp_numerical_stability() {
        std::array<size_t, 3> dims = {2, 2, 4};
        MultiVector<double, 3> zeros(dims);
        for (size_t i = 0; i < dims[0]; ++i)
            for (size_t j = 0; j < dims[1]; ++j)
                for (size_t k = 0; k < dims[2]; ++k)
                    zeros.at({i, j, k}) = 0.0;
        auto lse0 = zeros.logsumexp();
        ASSERT_LT(std::abs(lse0.at({0, 0}) - std::log(4.0)), 1e-10);

        MultiVector<double, 3> large(dims);
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 2; ++j)
                for (size_t k = 0; k < 4; ++k)
                    large.at({i, j, k}) = 700.0;
        auto lse_large = large.logsumexp();
        ASSERT_LT(std::abs(lse_large.at({0, 0}) - (700.0 + std::log(4.0))), 1e-6);
    }

    void test_softmax_normalization() {
        std::array<size_t, 3> dims = {2, 2, 5};
        MultiVector<double, 3> mv(dims);
        std::mt19937 gen(123);
        std::uniform_real_distribution<double> dis(-2.0, 2.0);
        for (size_t i = 0; i < dims[0]; ++i)
            for (size_t j = 0; j < dims[1]; ++j)
                for (size_t k = 0; k < dims[2]; ++k)
                    mv.at({i, j, k}) = dis(gen);
        auto sm = mv.softmax(false);
        for (size_t i = 0; i < dims[0]; ++i)
            for (size_t j = 0; j < dims[1]; ++j) {
                double row_sum = 0.0;
                for (size_t k = 0; k < dims[2]; ++k) row_sum += sm.at({i, j, k});
                ASSERT_LT(std::abs(row_sum - 1.0), 1e-10);
            }
    }

    void test_ragged_bounds_and_offsets() {
        std::array<size_t, 2> dims = {2, 2};
        std::vector<size_t> ragged = {3, 4};
        RaggedMultiVector<int, 3> rmv(dims, std::span<size_t const>(ragged));
        // For N=3: total_entries = product of first (N-2) dims = 2, data_size = sum(ragged) = 7
        ASSERT_EQ(rmv.total_size(), 2u * (3u + 4u));
        for (size_t i = 0; i < 2; ++i)
            for (size_t j = 0; j < 2; ++j) {
                size_t inner_len = (j == 0) ? 3u : 4u;
                size_t count = 0;
                for (auto it = rmv.inner_begin({i, j}); it != rmv.inner_end({i, j}); ++it) ++count;
                ASSERT_EQ(count, inner_len);
            }
    }

    void test_at_bounds_check() {
        std::array<size_t, 3> dims = {2, 3, 4};
        MultiVector<int, 3> mv(dims);
        ASSERT_THROWS(mv.at({2, 0, 0}));
        ASSERT_THROWS(mv.at({0, 4, 0}));
        ASSERT_THROWS(mv.at({0, 0, 5}));
        ASSERT_NO_THROW(mv.at({1, 2, 3}));
    }
};

} // namespace moire_test

int main() {
    auto suite = std::make_unique<moire_test::MultiVectorTestSuite>();
    moire_test::TestRegistry::register_suite(std::move(suite));
    moire_test::TestRegistry::run_all_tests();
    return 0;
}

