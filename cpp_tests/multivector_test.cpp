#include <iostream>
#include "multivector.h"
#include <vector>
#include <cassert>
#include <span>

void test_multivector_construction() {
    std::array<size_t, 3> dimensions = {2, 3, 4};
    MultiVector<int, 3> mv(dimensions);
    assert(mv.total_size() == 24);
    std::cout << "Construction test passed.\n";

    mv.inner_fill({0, 0}, 1);
    mv.inner_fill({0, 1}, 2);
    mv.inner_fill({0, 2}, 3);
    mv.inner_fill({1, 0}, 4);
    mv.inner_fill({1, 1}, 5);
    mv.inner_fill({1, 2}, 6);
    
    std::cout << "MultiVector: " << std::endl;
    for (size_t i = 0; i < mv.total_size(); ++i) {
        std::cout << mv.data()[i] << " ";
    }
    std::cout << std::endl;
}

void test_multivector_element_access() {
    std::array<size_t, 3> dimensions = {2, 3, 4};
    MultiVector<int, 3> mv(dimensions);
    mv.at({0, 0, 0}) = 42;
    mv.at({0, 0, 1}) = 43;
    mv.at({0, 0, 2}) = 44;
    mv.at({0, 0, 3}) = 45;
    mv.at({1, 2, 3}) = 46;

    std::cout << "MultiVector: " << std::endl;
    for (size_t i = 0; i < mv.total_size(); ++i) {
        std::cout << mv.data()[i] << " ";
    }
    std::cout << std::endl;
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

void test_multivector_reduce() {
    // Create a 2x3x4 MultiVector
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
    assert(sums.dimensions() == (std::array<size_t, 2>{2, 3}));
    assert(sums.at({0, 0}) == 10);  // 1+2+3+4
    assert(sums.at({0, 1}) == 14);  // 2+3+4+5
    assert(sums.at({0, 2}) == 18);  // 3+4+5+6
    assert(sums.at({1, 0}) == 22);  // 4+5+6+7
    assert(sums.at({1, 1}) == 26);  // 5+6+7+8
    assert(sums.at({1, 2}) == 30);  // 6+7+8+9

    // Test product reduction
    auto products = mv.product();
    assert(products.dimensions() == (std::array<size_t, 2>{2, 3}));
    assert(products.at({0, 0}) == 24);   // 1*2*3*4
    assert(products.at({0, 1}) == 120);  // 2*3*4*5
    assert(products.at({0, 2}) == 360);  // 3*4*5*6
    assert(products.at({1, 0}) == 840);  // 4*5*6*7
    assert(products.at({1, 1}) == 1680); // 5*6*7*8
    assert(products.at({1, 2}) == 3024); // 6*7*8*9

    // Test max reduction
    auto maxes = mv.max();
    assert(maxes.dimensions() == (std::array<size_t, 2>{2, 3}));
    assert(maxes.at({0, 0}) == 4);  // max(1,2,3,4)
    assert(maxes.at({0, 1}) == 5);  // max(2,3,4,5)
    assert(maxes.at({0, 2}) == 6);  // max(3,4,5,6)
    assert(maxes.at({1, 0}) == 7);  // max(4,5,6,7)
    assert(maxes.at({1, 1}) == 8);  // max(5,6,7,8)
    assert(maxes.at({1, 2}) == 9);  // max(6,7,8,9)

    // Test min reduction
    auto mins = mv.min();
    assert(mins.dimensions() == (std::array<size_t, 2>{2, 3}));
    assert(mins.at({0, 0}) == 1);  // min(1,2,3,4)
    assert(mins.at({0, 1}) == 2);  // min(2,3,4,5)
    assert(mins.at({0, 2}) == 3);  // min(3,4,5,6)
    assert(mins.at({1, 0}) == 4);  // min(4,5,6,7)
    assert(mins.at({1, 1}) == 5);  // min(5,6,7,8)
    assert(mins.at({1, 2}) == 6);  // min(6,7,8,9)

    std::cout << "Reduction test passed.\n";
}

// Test case for logsumexp reduction
void test_multivector_logsumexp() {
    std::array<size_t, 3> dimensions = {2, 3, 4};
    
    // Create a new MultiVector with floating-point values for logsumexp
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
    assert(logsumexps.dimensions() == (std::array<size_t, 2>{2, 3}));
    
    // Calculate expected values manually
    // For {0, 0}: log(exp(1) + exp(2) + exp(3) + exp(4)) = log(e + e^2 + e^3 + e^4) ≈ 4.4402
    assert(std::abs(logsumexps.at({0, 0}) - 4.4402) < 0.0001);
    
    // For {0, 1}: log(exp(2) + exp(3) + exp(4) + exp(5)) = log(e^2 + e^3 + e^4 + e^5) ≈ 5.4402
    assert(std::abs(logsumexps.at({0, 1}) - 5.4402) < 0.0001);
    
    // For {0, 2}: log(exp(3) + exp(4) + exp(5) + exp(6)) = log(e^3 + e^4 + e^5 + e^6) ≈ 6.4402
    assert(std::abs(logsumexps.at({0, 2}) - 6.4402) < 0.0001);
    
    // For {1, 0}: log(exp(4) + exp(5) + exp(6) + exp(7)) = log(e^4 + e^5 + e^6 + e^7) ≈ 7.4402
    assert(std::abs(logsumexps.at({1, 0}) - 7.4402) < 0.0001);
    
    // For {1, 1}: log(exp(5) + exp(6) + exp(7) + exp(8)) = log(e^5 + e^6 + e^7 + e^8) ≈ 8.4402
    assert(std::abs(logsumexps.at({1, 1}) - 8.4402) < 0.0001);
    
    // For {1, 2}: log(exp(6) + exp(7) + exp(8) + exp(9)) = log(e^6 + e^7 + e^8 + e^9) ≈ 9.4402
    assert(std::abs(logsumexps.at({1, 2}) - 9.4402) < 0.0001);
    
    // Test parallel_logsumexp
    auto parallel_logsumexps = mv_double.parallel_logsumexp();
    assert(parallel_logsumexps.dimensions() == (std::array<size_t, 2>{2, 3}));
    
    // Verify that parallel results match sequential results
    assert(std::abs(parallel_logsumexps.at({0, 0}) - logsumexps.at({0, 0})) < 0.0001);
    assert(std::abs(parallel_logsumexps.at({0, 1}) - logsumexps.at({0, 1})) < 0.0001);
    assert(std::abs(parallel_logsumexps.at({0, 2}) - logsumexps.at({0, 2})) < 0.0001);
    assert(std::abs(parallel_logsumexps.at({1, 0}) - logsumexps.at({1, 0})) < 0.0001);
    assert(std::abs(parallel_logsumexps.at({1, 1}) - logsumexps.at({1, 1})) < 0.0001);
    assert(std::abs(parallel_logsumexps.at({1, 2}) - logsumexps.at({1, 2})) < 0.0001);
    
    // Test with a more challenging case (values with large differences)
    MultiVector<double, 2> mv_challenge({1, 4});
    mv_challenge.at({0, 0}) = 100.0;
    mv_challenge.at({0, 1}) = 101.0;
    mv_challenge.at({0, 2}) = 102.0;
    mv_challenge.at({0, 3}) = 103.0;
    
    auto challenge_logsumexp = mv_challenge.logsumexp();
    // Expected: log(exp(100) + exp(101) + exp(102) + exp(103))
    // = 103 + log(e^-3 + e^-2 + e^-1 + 1)
    // ≈ 103 + log(1.553)
    // ≈ 103.44
    assert(std::abs(challenge_logsumexp.at({0}) - 103.44) < 0.001);
    
    std::cout << "Logsumexp reduction test passed.\n";
}

void test_multivector_transform() {
    std::cout << "Testing transform operations...\n";
    
    // Create a 2x3x4 MultiVector
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
    
    // Create a copy for testing
    MultiVector<int, 3> mv_copy = mv;
    
    // Create values for each inner slice (2*3 = 6 values)
    std::vector<int> add_values = {1, 2, 3, 4, 5, 6};
    
    // Test add operation
    mv = mv.add(add_values);
    
    // Verify results
    assert(mv.at({0, 0, 0}) == 2);  // 1+1
    assert(mv.at({0, 0, 1}) == 3);  // 2+1
    assert(mv.at({0, 0, 2}) == 4);  // 3+1
    assert(mv.at({0, 0, 3}) == 5);  // 4+1
    
    assert(mv.at({0, 1, 0}) == 4);  // 2+2
    assert(mv.at({0, 1, 1}) == 5);  // 3+2
    assert(mv.at({0, 1, 2}) == 6);  // 4+2
    assert(mv.at({0, 1, 3}) == 7);  // 5+2
    
    assert(mv.at({0, 2, 0}) == 6);  // 3+3
    assert(mv.at({0, 2, 1}) == 7);  // 4+3
    assert(mv.at({0, 2, 2}) == 8);  // 5+3
    assert(mv.at({0, 2, 3}) == 9);  // 6+3
    
    assert(mv.at({1, 0, 0}) == 8);  // 4+4
    assert(mv.at({1, 0, 1}) == 9);  // 5+4
    assert(mv.at({1, 0, 2}) == 10); // 6+4
    assert(mv.at({1, 0, 3}) == 11); // 7+4
    
    assert(mv.at({1, 1, 0}) == 10); // 5+5
    assert(mv.at({1, 1, 1}) == 11); // 6+5
    assert(mv.at({1, 1, 2}) == 12); // 7+5
    assert(mv.at({1, 1, 3}) == 13); // 8+5
    
    assert(mv.at({1, 2, 0}) == 12); // 6+6
    assert(mv.at({1, 2, 1}) == 13); // 7+6
    assert(mv.at({1, 2, 2}) == 14); // 8+6
    assert(mv.at({1, 2, 3}) == 15); // 9+6
    
    // Reset for next test
    mv = mv_copy;
    
    // Test subtract operation
    mv = mv.subtract(add_values);
    
    // Verify results
    assert(mv.at({0, 0, 0}) == 0);  // 1-1
    assert(mv.at({0, 0, 1}) == 1);  // 2-1
    assert(mv.at({0, 0, 2}) == 2);  // 3-1
    assert(mv.at({0, 0, 3}) == 3);  // 4-1
    
    assert(mv.at({0, 1, 0}) == 0);  // 2-2
    assert(mv.at({0, 1, 1}) == 1);  // 3-2
    assert(mv.at({0, 1, 2}) == 2);  // 4-2
    assert(mv.at({0, 1, 3}) == 3);  // 5-2
    
    assert(mv.at({0, 2, 0}) == 0);  // 3-3
    assert(mv.at({0, 2, 1}) == 1);  // 4-3
    assert(mv.at({0, 2, 2}) == 2);  // 5-3
    assert(mv.at({0, 2, 3}) == 3);  // 6-3
    
    assert(mv.at({1, 0, 0}) == 0);  // 4-4
    assert(mv.at({1, 0, 1}) == 1);  // 5-4
    assert(mv.at({1, 0, 2}) == 2);  // 6-4
    assert(mv.at({1, 0, 3}) == 3);  // 7-4
    
    assert(mv.at({1, 1, 0}) == 0);  // 5-5
    assert(mv.at({1, 1, 1}) == 1);  // 6-5
    assert(mv.at({1, 1, 2}) == 2);  // 7-5
    assert(mv.at({1, 1, 3}) == 3);  // 8-5
    
    assert(mv.at({1, 2, 0}) == 0);  // 6-6
    assert(mv.at({1, 2, 1}) == 1);  // 7-6
    assert(mv.at({1, 2, 2}) == 2);  // 8-6
    assert(mv.at({1, 2, 3}) == 3);  // 9-6
    
    // Reset for next test
    mv = mv_copy;
    
    // Test multiply operation
    std::vector<int> mult_values = {2, 3, 4, 5, 6, 7};
    mv = mv.multiply(mult_values);
    
    // Verify results
    assert(mv.at({0, 0, 0}) == 2);   // 1*2
    assert(mv.at({0, 0, 1}) == 4);   // 2*2
    assert(mv.at({0, 0, 2}) == 6);   // 3*2
    assert(mv.at({0, 0, 3}) == 8);   // 4*2
    
    assert(mv.at({0, 1, 0}) == 6);   // 2*3
    assert(mv.at({0, 1, 1}) == 9);   // 3*3
    assert(mv.at({0, 1, 2}) == 12);  // 4*3
    assert(mv.at({0, 1, 3}) == 15);  // 5*3
    
    assert(mv.at({0, 2, 0}) == 12);  // 3*4
    assert(mv.at({0, 2, 1}) == 16);  // 4*4
    assert(mv.at({0, 2, 2}) == 20);  // 5*4
    assert(mv.at({0, 2, 3}) == 24);  // 6*4
    
    assert(mv.at({1, 0, 0}) == 20);  // 4*5
    assert(mv.at({1, 0, 1}) == 25);  // 5*5
    assert(mv.at({1, 0, 2}) == 30);  // 6*5
    assert(mv.at({1, 0, 3}) == 35);  // 7*5
    
    assert(mv.at({1, 1, 0}) == 30);  // 5*6
    assert(mv.at({1, 1, 1}) == 36);  // 6*6
    assert(mv.at({1, 1, 2}) == 42);  // 7*6
    assert(mv.at({1, 1, 3}) == 48);  // 8*6
    
    assert(mv.at({1, 2, 0}) == 42);  // 6*7
    assert(mv.at({1, 2, 1}) == 49);  // 7*7
    assert(mv.at({1, 2, 2}) == 56);  // 8*7
    assert(mv.at({1, 2, 3}) == 63);  // 9*7
    
    // Reset for next test
    mv = mv_copy;
    
    // Test divide operation
    std::vector<int> div_values = {2, 2, 2, 2, 2, 2};
    mv = mv.divide(div_values);
    
    // Verify results (integer division)
    assert(mv.at({0, 0, 0}) == 0);  // 1/2
    assert(mv.at({0, 0, 1}) == 1);  // 2/2
    assert(mv.at({0, 0, 2}) == 1);  // 3/2
    assert(mv.at({0, 0, 3}) == 2);  // 4/2
    
    assert(mv.at({0, 1, 0}) == 1);  // 2/2
    assert(mv.at({0, 1, 1}) == 1);  // 3/2
    assert(mv.at({0, 1, 2}) == 2);  // 4/2
    assert(mv.at({0, 1, 3}) == 2);  // 5/2
    
    assert(mv.at({0, 2, 0}) == 1);  // 3/2
    assert(mv.at({0, 2, 1}) == 2);  // 4/2
    assert(mv.at({0, 2, 2}) == 2);  // 5/2
    assert(mv.at({0, 2, 3}) == 3);  // 6/2
    
    assert(mv.at({1, 0, 0}) == 2);  // 4/2
    assert(mv.at({1, 0, 1}) == 2);  // 5/2
    assert(mv.at({1, 0, 2}) == 3);  // 6/2
    assert(mv.at({1, 0, 3}) == 3);  // 7/2
    
    assert(mv.at({1, 1, 0}) == 2);  // 5/2
    assert(mv.at({1, 1, 1}) == 3);  // 6/2
    assert(mv.at({1, 1, 2}) == 3);  // 7/2
    assert(mv.at({1, 1, 3}) == 4);  // 8/2
    
    assert(mv.at({1, 2, 0}) == 3);  // 6/2
    assert(mv.at({1, 2, 1}) == 3);  // 7/2
    assert(mv.at({1, 2, 2}) == 4);  // 8/2
    assert(mv.at({1, 2, 3}) == 4);  // 9/2
    
    // Test parallel versions using std::execution::par
    // Reset for parallel tests
    mv = mv_copy;
    
    // Test parallel_add
    mv = mv.add(add_values, std::execution::par);
    
    // Verify results (same as sequential add)
    assert(mv.at({0, 0, 0}) == 2);  // 1+1
    assert(mv.at({0, 0, 1}) == 3);  // 2+1
    assert(mv.at({0, 0, 2}) == 4);  // 3+1
    assert(mv.at({0, 0, 3}) == 5);  // 4+1
    
    // Reset for next parallel test
    mv = mv_copy;
    
    // Test parallel_multiply
    mv = mv.multiply(mult_values, std::execution::par);
    
    // Verify results (same as sequential multiply)
    assert(mv.at({0, 0, 0}) == 2);   // 1*2
    assert(mv.at({0, 0, 1}) == 4);   // 2*2
    assert(mv.at({0, 0, 2}) == 6);   // 3*2
    assert(mv.at({0, 0, 3}) == 8);   // 4*2
    
    std::cout << "Transform operations test passed.\n";
}

void test_multivector_element_transform() {
    std::cout << "Testing element_transform operations...\n";
    
    // Create a 2x3x4 MultiVector
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
    
    // Create a copy for testing
    MultiVector<int, 3> mv_copy = mv;
    
    // Create values for each element position in the inner dimension (4 values)
    std::vector<int> add_values = {1, 2, 3, 4};
    
    // Test element_add operation
    mv = mv.element_add(add_values);
    
    // Verify results
    // First slice (i=0)
    assert(mv.at({0, 0, 0}) == 2);  // 1+1
    assert(mv.at({0, 0, 1}) == 4);  // 2+2
    assert(mv.at({0, 0, 2}) == 6);  // 3+3
    assert(mv.at({0, 0, 3}) == 8);  // 4+4
    
    assert(mv.at({0, 1, 0}) == 3);  // 2+1
    assert(mv.at({0, 1, 1}) == 5);  // 3+2
    assert(mv.at({0, 1, 2}) == 7);  // 4+3
    assert(mv.at({0, 1, 3}) == 9);  // 5+4
    
    assert(mv.at({0, 2, 0}) == 4);  // 3+1
    assert(mv.at({0, 2, 1}) == 6);  // 4+2
    assert(mv.at({0, 2, 2}) == 8);  // 5+3
    assert(mv.at({0, 2, 3}) == 10); // 6+4
    
    // Second slice (i=1)
    assert(mv.at({1, 0, 0}) == 5);  // 4+1
    assert(mv.at({1, 0, 1}) == 7);  // 5+2
    assert(mv.at({1, 0, 2}) == 9);  // 6+3
    assert(mv.at({1, 0, 3}) == 11); // 7+4
    
    assert(mv.at({1, 1, 0}) == 6);  // 5+1
    assert(mv.at({1, 1, 1}) == 8);  // 6+2
    assert(mv.at({1, 1, 2}) == 10); // 7+3
    assert(mv.at({1, 1, 3}) == 12); // 8+4
    
    assert(mv.at({1, 2, 0}) == 7);  // 6+1
    assert(mv.at({1, 2, 1}) == 9);  // 7+2
    assert(mv.at({1, 2, 2}) == 11); // 8+3
    assert(mv.at({1, 2, 3}) == 13); // 9+4
    
    // Reset for next test
    mv = mv_copy;
    
    // Test element_subtract operation
    mv = mv.element_subtract(add_values);
    
    // Verify results
    // First slice (i=0)
    assert(mv.at({0, 0, 0}) == 0);  // 1-1
    assert(mv.at({0, 0, 1}) == 0);  // 2-2
    assert(mv.at({0, 0, 2}) == 0);  // 3-3
    assert(mv.at({0, 0, 3}) == 0);  // 4-4
    
    assert(mv.at({0, 1, 0}) == 1);  // 2-1
    assert(mv.at({0, 1, 1}) == 1);  // 3-2
    assert(mv.at({0, 1, 2}) == 1);  // 4-3
    assert(mv.at({0, 1, 3}) == 1);  // 5-4
    
    assert(mv.at({0, 2, 0}) == 2);  // 3-1
    assert(mv.at({0, 2, 1}) == 2);  // 4-2
    assert(mv.at({0, 2, 2}) == 2);  // 5-3
    assert(mv.at({0, 2, 3}) == 2);  // 6-4
    
    // Second slice (i=1)
    assert(mv.at({1, 0, 0}) == 3);  // 4-1
    assert(mv.at({1, 0, 1}) == 3);  // 5-2
    assert(mv.at({1, 0, 2}) == 3);  // 6-3
    assert(mv.at({1, 0, 3}) == 3);  // 7-4
    
    assert(mv.at({1, 1, 0}) == 4);  // 5-1
    assert(mv.at({1, 1, 1}) == 4);  // 6-2
    assert(mv.at({1, 1, 2}) == 4);  // 7-3
    assert(mv.at({1, 1, 3}) == 4);  // 8-4
    
    assert(mv.at({1, 2, 0}) == 5);  // 6-1
    assert(mv.at({1, 2, 1}) == 5);  // 7-2
    assert(mv.at({1, 2, 2}) == 5);  // 8-3
    assert(mv.at({1, 2, 3}) == 5);  // 9-4
    
    // Reset for next test
    mv = mv_copy;
    
    // Test element_multiply operation
    std::vector<int> mult_values = {2, 3, 4, 5};
    mv = mv.element_multiply(mult_values);
    
    // Verify results
    // First slice (i=0)
    assert(mv.at({0, 0, 0}) == 2);   // 1*2
    assert(mv.at({0, 0, 1}) == 6);   // 2*3
    assert(mv.at({0, 0, 2}) == 12);  // 3*4
    assert(mv.at({0, 0, 3}) == 20);  // 4*5
    
    assert(mv.at({0, 1, 0}) == 4);   // 2*2
    assert(mv.at({0, 1, 1}) == 9);   // 3*3
    assert(mv.at({0, 1, 2}) == 16);  // 4*4
    assert(mv.at({0, 1, 3}) == 25);  // 5*5
    
    assert(mv.at({0, 2, 0}) == 6);   // 3*2
    assert(mv.at({0, 2, 1}) == 12);  // 4*3
    assert(mv.at({0, 2, 2}) == 20);  // 5*4
    assert(mv.at({0, 2, 3}) == 30);  // 6*5
    
    // Second slice (i=1)
    assert(mv.at({1, 0, 0}) == 8);   // 4*2
    assert(mv.at({1, 0, 1}) == 15);  // 5*3
    assert(mv.at({1, 0, 2}) == 24);  // 6*4
    assert(mv.at({1, 0, 3}) == 35);  // 7*5
    
    assert(mv.at({1, 1, 0}) == 10);  // 5*2
    assert(mv.at({1, 1, 1}) == 18);  // 6*3
    assert(mv.at({1, 1, 2}) == 28);  // 7*4
    assert(mv.at({1, 1, 3}) == 40);  // 8*5
    
    assert(mv.at({1, 2, 0}) == 12);  // 6*2
    assert(mv.at({1, 2, 1}) == 21);  // 7*3
    assert(mv.at({1, 2, 2}) == 32);  // 8*4
    assert(mv.at({1, 2, 3}) == 45);  // 9*5
    
    // Reset for next test
    mv = mv_copy;
    
    // Test element_divide operation
    std::vector<int> div_values = {2, 2, 2, 2};
    mv = mv.element_divide(div_values);
    
    // Verify results (integer division)
    // First slice (i=0)
    assert(mv.at({0, 0, 0}) == 0);  // 1/2
    assert(mv.at({0, 0, 1}) == 1);  // 2/2
    assert(mv.at({0, 0, 2}) == 1);  // 3/2
    assert(mv.at({0, 0, 3}) == 2);  // 4/2
    
    assert(mv.at({0, 1, 0}) == 1);  // 2/2
    assert(mv.at({0, 1, 1}) == 1);  // 3/2
    assert(mv.at({0, 1, 2}) == 2);  // 4/2
    assert(mv.at({0, 1, 3}) == 2);  // 5/2
    
    assert(mv.at({0, 2, 0}) == 1);  // 3/2
    assert(mv.at({0, 2, 1}) == 2);  // 4/2
    assert(mv.at({0, 2, 2}) == 2);  // 5/2
    assert(mv.at({0, 2, 3}) == 3);  // 6/2
    
    // Second slice (i=1)
    assert(mv.at({1, 0, 0}) == 2);  // 4/2
    assert(mv.at({1, 0, 1}) == 2);  // 5/2
    assert(mv.at({1, 0, 2}) == 3);  // 6/2
    assert(mv.at({1, 0, 3}) == 3);  // 7/2
    
    assert(mv.at({1, 1, 0}) == 2);  // 5/2
    assert(mv.at({1, 1, 1}) == 3);  // 6/2
    assert(mv.at({1, 1, 2}) == 3);  // 7/2
    assert(mv.at({1, 1, 3}) == 4);  // 8/2
    
    assert(mv.at({1, 2, 0}) == 3);  // 6/2
    assert(mv.at({1, 2, 1}) == 3);  // 7/2
    assert(mv.at({1, 2, 2}) == 4);  // 8/2
    assert(mv.at({1, 2, 3}) == 4);  // 9/2
    
    // Test parallel versions
    // Reset for parallel tests
    mv = mv_copy;
    
    // Test parallel_element_add
    mv = mv.parallel_element_add(add_values);
    
    // Verify results (same as sequential add)
    assert(mv.at({0, 0, 0}) == 2);  // 1+1
    assert(mv.at({0, 0, 1}) == 4);  // 2+2
    assert(mv.at({0, 0, 2}) == 6);  // 3+3
    assert(mv.at({0, 0, 3}) == 8);  // 4+4
    
    // Reset for next parallel test
    mv = mv_copy;
    
    // Test parallel_element_multiply
    mv = mv.parallel_element_multiply(mult_values);
    
    // Verify results (same as sequential multiply)
    assert(mv.at({0, 0, 0}) == 2);   // 1*2
    assert(mv.at({0, 0, 1}) == 6);   // 2*3
    assert(mv.at({0, 0, 2}) == 12);  // 3*4
    assert(mv.at({0, 0, 3}) == 20);  // 4*5
    
    // Test N=1 specialization
    MultiVector<int, 1> mv1d({4});
    mv1d.at({0}) = 1;
    mv1d.at({1}) = 2;
    mv1d.at({2}) = 3;
    mv1d.at({3}) = 4;
    
    std::vector<int> values1d = {2, 3, 4, 5};
    
    // Test element_add
    mv1d = mv1d.element_add(values1d);
    assert(mv1d.at({0}) == 3);  // 1+2
    assert(mv1d.at({1}) == 5);  // 2+3
    assert(mv1d.at({2}) == 7);  // 3+4
    assert(mv1d.at({3}) == 9);  // 4+5
    
    // Test parallel_element_add
    MultiVector<int, 1> mv1d_copy = mv1d;
    mv1d_copy = mv1d_copy.parallel_element_add(values1d);
    assert(mv1d_copy.at({0}) == 5);  // 3+2
    assert(mv1d_copy.at({1}) == 8);  // 5+3
    assert(mv1d_copy.at({2}) == 11); // 7+4
    assert(mv1d_copy.at({3}) == 14); // 9+5
    
    std::cout << "Element transform operations test passed.\n";
}

// Test case for RaggedMultiVector construction
void test_ragged_multivector_construction() {
    std::array<size_t, 3> dimensions = {2, 3, 3};
    std::vector<size_t> ragged_dimensions = {2, 3, 4};
    RaggedMultiVector<int, 4> raggedMultiVector(dimensions, std::span<size_t const>(ragged_dimensions));
    // RaggedMultiVector<int, 1> raggedMultiVector2(std::array<size_t, 0>{}, std::span<size_t const>(ragged_dimensions));
    
    // Check if the total size is correct
    size_t expected_total_size = 2 * 3 * (2 + 3 + 4);
    assert(raggedMultiVector.total_size() == expected_total_size);
    // assert(raggedMultiVector2.total_size() == 9);

    // Check dimensions
    auto dims = raggedMultiVector.dimensions();
    assert(dims[0] == 2);
    assert(dims[1] == 3);
    assert(dims[2] == 3); // The last dimension is the size of ragged_dimensions
    
    // auto dims2 = raggedMultiVector2.dimensions();
    // assert(dims2[0] == 3); // The size of ragged_dimensions

    std::cout << "RaggedMultiVector construction test passed." << std::endl;
}

// Test case for RaggedMultiVector element access
void test_ragged_multivector_element_access() {
    std::array<size_t, 2> dimensions = {2, 3}; // 2 samples, 3 loci, for example
    std::vector<size_t> ragged_dimensions = {2, 3, 4};
    RaggedMultiVector<int, 3> raggedMultiVector(dimensions, std::span<size_t const>(ragged_dimensions));
    // RaggedMultiVector<float, 1> raggedMultiVector2(std::array<size_t, 0>{}, std::span<size_t const>(ragged_dimensions));

    // Access elements using the at method
    raggedMultiVector.at({0, 0, 0}) = 10;
    raggedMultiVector.at({0, 1, 1}) = 20;
    raggedMultiVector.at({1, 2, 2}) = 30;
    raggedMultiVector.at({0, 2, 0}) = 40;
    raggedMultiVector.at({1, 0, 1}) = 50;
    raggedMultiVector.at({1, 1, 2}) = 60;

    // raggedMultiVector2.at({0}) = 10.0f;
    // raggedMultiVector2.at({1}) = 20.0f;
    // raggedMultiVector2.at({2}) = 30.0f; 

    // Verify the values
    assert(raggedMultiVector.at({0, 0, 0}) == 10);
    assert(raggedMultiVector.at({0, 1, 1}) == 20);
    assert(raggedMultiVector.at({1, 2, 2}) == 30);
    assert(raggedMultiVector.at({0, 2, 0}) == 40);
    assert(raggedMultiVector.at({1, 0, 1}) == 50);
    assert(raggedMultiVector.at({1, 1, 2}) == 60);

    // assert(raggedMultiVector2.at({0}) == 10.0f);
    // assert(raggedMultiVector2.at({1}) == 20.0f);
    // assert(raggedMultiVector2.at({2}) == 30.0f);

    // Check total size
    assert(raggedMultiVector.total_size() == 18);
    // assert(raggedMultiVector2.total_size() == 9);

    std::cout << "RaggedMultiVector element access test passed." << std::endl;
}

// Test case for RaggedMultiVector iteration
void test_ragged_multivector_iteration() {
    std::array<size_t, 2> dimensions = {2, 3};
    std::vector<size_t> ragged_dimensions = {2, 3, 4};
    RaggedMultiVector<int, 3> raggedMultiVector(dimensions, std::span<size_t const>(ragged_dimensions));

    // Fill the first inner dimension with values
    for (size_t i = 0; i < ragged_dimensions[0]; ++i) {
        raggedMultiVector.at({0, 0, i}) = i;
    }

    // Iterate over the inner dimension and verify values
    int index = 0;
    for (auto it = raggedMultiVector.inner_begin({0, 0}); it != raggedMultiVector.inner_end({0, 0}); ++it) {
        assert(*it == index++);
    }

    std::cout << "RaggedMultiVector iteration test passed." << std::endl;
}

void test_multivector_n1_specialization() {
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
    assert(sum == 15);  // 1+2+3+4+5
    
    // Test product reduction - should return a scalar
    int product = mv.product();
    assert(product == 120);  // 1*2*3*4*5
    
    // Test max reduction - should return a scalar
    int max_val = mv.max();
    assert(max_val == 5);  // max(1,2,3,4,5)
    
    // Test min reduction - should return a scalar
    int min_val = mv.min();
    assert(min_val == 1);  // min(1,2,3,4,5)
    
    // Test custom reduction - should return a scalar
    auto custom_reduce = mv.reduce([](int a, int b) { 
        return a + b * 2; 
    }, 0);
    assert(custom_reduce == 30);  // 0 + 1*2 + 2*2 + 3*2 + 4*2 + 5*2 = 0 + 2 + 4 + 6 + 8 + 10 = 30
    
    // Test parallel reductions
    int parallel_sum = mv.parallel_sum();
    assert(parallel_sum == 15);
    
    int parallel_product = mv.parallel_product();
    assert(parallel_product == 120);
    
    int parallel_max = mv.parallel_max();
    assert(parallel_max == 5);
    
    int parallel_min = mv.parallel_min();
    assert(parallel_min == 1);
    
    // Test logsumexp for floating point
    MultiVector<double, 1> mv_double({5});
    mv_double.at({0}) = 1.0;
    mv_double.at({1}) = 2.0;
    mv_double.at({2}) = 3.0;
    mv_double.at({3}) = 4.0;
    mv_double.at({4}) = 5.0;
    
    double logsumexp_val = mv_double.logsumexp();
    // log(sum(exp(x))) = log(exp(1) + exp(2) + exp(3) + exp(4) + exp(5))
    // This is approximately 5.451914
    assert(std::abs(logsumexp_val - 5.451914) < 0.0001);
    
    double parallel_logsumexp_val = mv_double.parallel_logsumexp();
    assert(std::abs(parallel_logsumexp_val - 5.451914) < 0.0001);
    
    std::cout << "N=1 MultiVector specialization test passed." << std::endl;
}

void test_multivector_elementwise_operations() {
    std::cout << "Testing elementwise operations between multivectors...\n";
    
    // Create two 2x3x4 MultiVectors
    std::array<size_t, 3> dimensions = {2, 3, 4};
    MultiVector<int, 3> mv1(dimensions);
    MultiVector<int, 3> mv2(dimensions);
    
    // Fill with test data
    // First slice (i=0)
    mv1.at({0, 0, 0}) = 1;  mv1.at({0, 0, 1}) = 2;  mv1.at({0, 0, 2}) = 3;  mv1.at({0, 0, 3}) = 4;
    mv1.at({0, 1, 0}) = 2;  mv1.at({0, 1, 1}) = 3;  mv1.at({0, 1, 2}) = 4;  mv1.at({0, 1, 3}) = 5;
    mv1.at({0, 2, 0}) = 3;  mv1.at({0, 2, 1}) = 4;  mv1.at({0, 2, 2}) = 5;  mv1.at({0, 2, 3}) = 6;
    
    // Second slice (i=1)
    mv1.at({1, 0, 0}) = 4;  mv1.at({1, 0, 1}) = 5;  mv1.at({1, 0, 2}) = 6;  mv1.at({1, 0, 3}) = 7;
    mv1.at({1, 1, 0}) = 5;  mv1.at({1, 1, 1}) = 6;  mv1.at({1, 1, 2}) = 7;  mv1.at({1, 1, 3}) = 8;
    mv1.at({1, 2, 0}) = 6;  mv1.at({1, 2, 1}) = 7;  mv1.at({1, 2, 2}) = 8;  mv1.at({1, 2, 3}) = 9;
    
    // Fill second multivector with different values
    mv2.at({0, 0, 0}) = 10;  mv2.at({0, 0, 1}) = 11;  mv2.at({0, 0, 2}) = 12;  mv2.at({0, 0, 3}) = 13;
    mv2.at({0, 1, 0}) = 11;  mv2.at({0, 1, 1}) = 12;  mv2.at({0, 1, 2}) = 13;  mv2.at({0, 1, 3}) = 14;
    mv2.at({0, 2, 0}) = 12;  mv2.at({0, 2, 1}) = 13;  mv2.at({0, 2, 2}) = 14;  mv2.at({0, 2, 3}) = 15;
    
    // Second slice (i=1)
    mv2.at({1, 0, 0}) = 13;  mv2.at({1, 0, 1}) = 14;  mv2.at({1, 0, 2}) = 15;  mv2.at({1, 0, 3}) = 16;
    mv2.at({1, 1, 0}) = 14;  mv2.at({1, 1, 1}) = 15;  mv2.at({1, 1, 2}) = 16;  mv2.at({1, 1, 3}) = 17;
    mv2.at({1, 2, 0}) = 15;  mv2.at({1, 2, 1}) = 16;  mv2.at({1, 2, 2}) = 17;  mv2.at({1, 2, 3}) = 18;
    
    // Test addition
    MultiVector<int, 3> sum = mv1 + mv2;
    
    // Verify results
    assert(sum.at({0, 0, 0}) == 11);  // 1+10
    assert(sum.at({0, 0, 1}) == 13);  // 2+11
    assert(sum.at({0, 0, 2}) == 15);  // 3+12
    assert(sum.at({0, 0, 3}) == 17);  // 4+13
    
    assert(sum.at({0, 1, 0}) == 13);  // 2+11
    assert(sum.at({0, 1, 1}) == 15);  // 3+12
    assert(sum.at({0, 1, 2}) == 17);  // 4+13
    assert(sum.at({0, 1, 3}) == 19);  // 5+14
    
    assert(sum.at({0, 2, 0}) == 15);  // 3+12
    assert(sum.at({0, 2, 1}) == 17);  // 4+13
    assert(sum.at({0, 2, 2}) == 19);  // 5+14
    assert(sum.at({0, 2, 3}) == 21);  // 6+15
    
    assert(sum.at({1, 0, 0}) == 17);  // 4+13
    assert(sum.at({1, 0, 1}) == 19);  // 5+14
    assert(sum.at({1, 0, 2}) == 21);  // 6+15
    assert(sum.at({1, 0, 3}) == 23);  // 7+16
    
    assert(sum.at({1, 1, 0}) == 19);  // 5+14
    assert(sum.at({1, 1, 1}) == 21);  // 6+15
    assert(sum.at({1, 1, 2}) == 23);  // 7+16
    assert(sum.at({1, 1, 3}) == 25);  // 8+17
    
    assert(sum.at({1, 2, 0}) == 21);  // 6+15
    assert(sum.at({1, 2, 1}) == 23);  // 7+16
    assert(sum.at({1, 2, 2}) == 25);  // 8+17
    assert(sum.at({1, 2, 3}) == 27);  // 9+18
    
    // Test subtraction
    MultiVector<int, 3> diff = mv2 - mv1;
    
    // Verify results
    assert(diff.at({0, 0, 0}) == 9);   // 10-1
    assert(diff.at({0, 0, 1}) == 9);   // 11-2
    assert(diff.at({0, 0, 2}) == 9);   // 12-3
    assert(diff.at({0, 0, 3}) == 9);   // 13-4
    
    assert(diff.at({0, 1, 0}) == 9);   // 11-2
    assert(diff.at({0, 1, 1}) == 9);   // 12-3
    assert(diff.at({0, 1, 2}) == 9);   // 13-4
    assert(diff.at({0, 1, 3}) == 9);   // 14-5
    
    assert(diff.at({0, 2, 0}) == 9);   // 12-3
    assert(diff.at({0, 2, 1}) == 9);   // 13-4
    assert(diff.at({0, 2, 2}) == 9);   // 14-5
    assert(diff.at({0, 2, 3}) == 9);   // 15-6
    
    assert(diff.at({1, 0, 0}) == 9);   // 13-4
    assert(diff.at({1, 0, 1}) == 9);   // 14-5
    assert(diff.at({1, 0, 2}) == 9);   // 15-6
    assert(diff.at({1, 0, 3}) == 9);   // 16-7
    
    assert(diff.at({1, 1, 0}) == 9);   // 14-5
    assert(diff.at({1, 1, 1}) == 9);   // 15-6
    assert(diff.at({1, 1, 2}) == 9);   // 16-7
    assert(diff.at({1, 1, 3}) == 9);   // 17-8
    
    assert(diff.at({1, 2, 0}) == 9);   // 15-6
    assert(diff.at({1, 2, 1}) == 9);   // 16-7
    assert(diff.at({1, 2, 2}) == 9);   // 17-8
    assert(diff.at({1, 2, 3}) == 9);   // 18-9
    
    // Test multiplication
    MultiVector<int, 3> product = mv1 * mv2;
    
    // Verify results
    assert(product.at({0, 0, 0}) == 10);   // 1*10
    assert(product.at({0, 0, 1}) == 22);   // 2*11
    assert(product.at({0, 0, 2}) == 36);   // 3*12
    assert(product.at({0, 0, 3}) == 52);   // 4*13
    
    assert(product.at({0, 1, 0}) == 22);   // 2*11
    assert(product.at({0, 1, 1}) == 36);   // 3*12
    assert(product.at({0, 1, 2}) == 52);   // 4*13
    assert(product.at({0, 1, 3}) == 70);   // 5*14
    
    assert(product.at({0, 2, 0}) == 36);   // 3*12
    assert(product.at({0, 2, 1}) == 52);   // 4*13
    assert(product.at({0, 2, 2}) == 70);   // 5*14
    assert(product.at({0, 2, 3}) == 90);   // 6*15
    
    assert(product.at({1, 0, 0}) == 52);   // 4*13
    assert(product.at({1, 0, 1}) == 70);   // 5*14
    assert(product.at({1, 0, 2}) == 90);   // 6*15
    assert(product.at({1, 0, 3}) == 112);  // 7*16
    
    assert(product.at({1, 1, 0}) == 70);   // 5*14
    assert(product.at({1, 1, 1}) == 90);   // 6*15
    assert(product.at({1, 1, 2}) == 112);  // 7*16
    assert(product.at({1, 1, 3}) == 136);  // 8*17
    
    assert(product.at({1, 2, 0}) == 90);   // 6*15
    assert(product.at({1, 2, 1}) == 112);  // 7*16
    assert(product.at({1, 2, 2}) == 136);  // 8*17
    assert(product.at({1, 2, 3}) == 162);  // 9*18
    
    // Test division (integer division)
    MultiVector<int, 3> quotient = mv2 / mv1;
    
    // Verify results
    assert(quotient.at({0, 0, 0}) == 10);  // 10/1
    assert(quotient.at({0, 0, 1}) == 5);   // 11/2
    assert(quotient.at({0, 0, 2}) == 4);   // 12/3
    assert(quotient.at({0, 0, 3}) == 3);   // 13/4
    
    assert(quotient.at({0, 1, 0}) == 5);   // 11/2
    assert(quotient.at({0, 1, 1}) == 4);   // 12/3
    assert(quotient.at({0, 1, 2}) == 3);   // 13/4
    assert(quotient.at({0, 1, 3}) == 2);   // 14/5
    
    assert(quotient.at({0, 2, 0}) == 4);   // 12/3
    assert(quotient.at({0, 2, 1}) == 3);   // 13/4
    assert(quotient.at({0, 2, 2}) == 2);   // 14/5
    assert(quotient.at({0, 2, 3}) == 2);   // 15/6
    
    assert(quotient.at({1, 0, 0}) == 3);   // 13/4
    assert(quotient.at({1, 0, 1}) == 2);   // 14/5
    assert(quotient.at({1, 0, 2}) == 2);   // 15/6
    assert(quotient.at({1, 0, 3}) == 2);   // 16/7
    
    assert(quotient.at({1, 1, 0}) == 2);   // 14/5
    assert(quotient.at({1, 1, 1}) == 2);   // 15/6
    assert(quotient.at({1, 1, 2}) == 2);   // 16/7
    assert(quotient.at({1, 1, 3}) == 2);   // 17/8
    
    assert(quotient.at({1, 2, 0}) == 2);   // 15/6
    assert(quotient.at({1, 2, 1}) == 2);   // 16/7
    assert(quotient.at({1, 2, 2}) == 2);   // 17/8
    assert(quotient.at({1, 2, 3}) == 2);   // 18/9
    
    // Test in-place operations
    MultiVector<int, 3> mv3 = mv1;  // Make a copy
    
    // Test +=
    mv3 += mv2;
    
    // Verify results (should be the same as sum)
    assert(mv3.at({0, 0, 0}) == 11);  // 1+10
    assert(mv3.at({0, 0, 1}) == 13);  // 2+11
    assert(mv3.at({0, 0, 2}) == 15);  // 3+12
    assert(mv3.at({0, 0, 3}) == 17);  // 4+13
    
    // Test -=
    mv3 = mv1;  // Reset
    mv3 -= mv2;
    
    // Verify results (should be the same as diff but with opposite sign)
    assert(mv3.at({0, 0, 0}) == -9);  // 1-10
    assert(mv3.at({0, 0, 1}) == -9);  // 2-11
    assert(mv3.at({0, 0, 2}) == -9);  // 3-12
    assert(mv3.at({0, 0, 3}) == -9);  // 4-13
    
    // Test *=
    mv3 = mv1;  // Reset
    mv3 *= mv2;
    
    // Verify results (should be the same as product)
    assert(mv3.at({0, 0, 0}) == 10);   // 1*10
    assert(mv3.at({0, 0, 1}) == 22);   // 2*11
    assert(mv3.at({0, 0, 2}) == 36);   // 3*12
    assert(mv3.at({0, 0, 3}) == 52);   // 4*13
    
    // Test /=
    mv3 = mv1;  // Reset
    mv3 /= mv2;
    
    // Verify results (should be the same as quotient but with opposite numerator/denominator)
    assert(mv3.at({0, 0, 0}) == 0);   // 1/10 (integer division)
    assert(mv3.at({0, 0, 1}) == 0);   // 2/11
    assert(mv3.at({0, 0, 2}) == 0);   // 3/12
    assert(mv3.at({0, 0, 3}) == 0);   // 4/13
    
    // Test with different dimensions (should throw)
    std::array<size_t, 3> different_dimensions = {2, 4, 4};
    MultiVector<int, 3> mv4(different_dimensions);
    
    bool exception_thrown = false;
    try {
        mv1 + mv4;  // This should throw
    } catch (const std::invalid_argument& e) {
        exception_thrown = true;
    }
    assert(exception_thrown);
    
    std::cout << "Elementwise operations test passed.\n";
}

void test_multivector_softmax() {
    std::cout << "Testing softmax operations...\n";
    
    // Create a 2x3x4 MultiVector with floating point values
    std::array<size_t, 3> dimensions = {2, 3, 4};
    MultiVector<double, 3> mv(dimensions);
    
    // Fill with test data
    // First slice (i=0)
    mv.at({0, 0, 0}) = 1.0;  mv.at({0, 0, 1}) = 2.0;  mv.at({0, 0, 2}) = 3.0;  mv.at({0, 0, 3}) = 4.0;
    mv.at({0, 1, 0}) = 2.0;  mv.at({0, 1, 1}) = 3.0;  mv.at({0, 1, 2}) = 4.0;  mv.at({0, 1, 3}) = 5.0;
    mv.at({0, 2, 0}) = 3.0;  mv.at({0, 2, 1}) = 4.0;  mv.at({0, 2, 2}) = 5.0;  mv.at({0, 2, 3}) = 6.0;
    
    // Second slice (i=1)
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
    assert(std::abs(sum0 - 1.0) < 1e-10);  // Sum should be 1.0
    
    // Verify results for second inner slice
    double sum1 = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        sum1 += softmax_result.at({0, 1, k});
    }
    assert(std::abs(sum1 - 1.0) < 1e-10);  // Sum should be 1.0
    
    // Test log-space softmax
    auto log_softmax_result = mv.softmax(true);
    
    // Verify results for first inner slice
    double log_sum0 = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        log_sum0 += std::exp(log_softmax_result.at({0, 0, k}));
    }
    assert(std::abs(log_sum0 - 1.0) < 1e-10);  // Sum should be 1.0
    
    // Verify results for second inner slice
    double log_sum1 = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        log_sum1 += std::exp(log_softmax_result.at({0, 1, k}));
    }
    assert(std::abs(log_sum1 - 1.0) < 1e-10);  // Sum should be 1.0
    
    // Test parallel softmax
    auto parallel_softmax_result = mv.parallel_softmax(false);
    
    // Verify results for first inner slice
    double parallel_sum0 = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        parallel_sum0 += parallel_softmax_result.at({0, 0, k});
    }
    assert(std::abs(parallel_sum0 - 1.0) < 1e-10);  // Sum should be 1.0
    
    // Verify results for second inner slice
    double parallel_sum1 = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        parallel_sum1 += parallel_softmax_result.at({0, 1, k});
    }
    assert(std::abs(parallel_sum1 - 1.0) < 1e-10);  // Sum should be 1.0
    
    // Test parallel log-space softmax
    auto parallel_log_softmax_result = mv.parallel_softmax(true);
    
    // Verify results for first inner slice
    double parallel_log_sum0 = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        parallel_log_sum0 += std::exp(parallel_log_softmax_result.at({0, 0, k}));
    }
    assert(std::abs(parallel_log_sum0 - 1.0) < 1e-10);  // Sum should be 1.0
    
    // Verify results for second inner slice
    double parallel_log_sum1 = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        parallel_log_sum1 += std::exp(parallel_log_softmax_result.at({0, 1, k}));
    }
    assert(std::abs(parallel_log_sum1 - 1.0) < 1e-10);  // Sum should be 1.0
    
    // Test with a more challenging case (values with large differences)
    MultiVector<double, 2> mv_challenge({1, 4});
    mv_challenge.at({0, 0}) = 100.0;
    mv_challenge.at({0, 1}) = 101.0;
    mv_challenge.at({0, 2}) = 102.0;
    mv_challenge.at({0, 3}) = 103.0;
    
    auto challenge_softmax = mv_challenge.softmax(false);
    auto challenge_log_softmax = mv_challenge.softmax(true);
    
    // Verify that the results are numerically stable
    double challenge_sum = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        challenge_sum += challenge_softmax.at({0, k});
    }
    assert(std::abs(challenge_sum - 1.0) < 1e-10);
    
    double challenge_log_sum = 0.0;
    for (size_t k = 0; k < 4; ++k) {
        challenge_log_sum += std::exp(challenge_log_softmax.at({0, k}));
    }
    assert(std::abs(challenge_log_sum - 1.0) < 1e-10);
    
    std::cout << "Softmax operations test passed.\n";
}

int main() {
    test_multivector_construction();
    test_multivector_element_access();
    test_multivector_dimension_size();
    test_multivector_iteration();
    test_multivector_reduce();
    test_multivector_transform();
    test_multivector_element_transform();
    test_multivector_logsumexp();
    test_ragged_multivector_construction();
    test_ragged_multivector_element_access();
    test_ragged_multivector_iteration();
    test_multivector_n1_specialization();
    test_multivector_elementwise_operations();
    std::cout << "All tests passed." << std::endl;
    return 0;
}


