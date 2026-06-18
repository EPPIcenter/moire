#include "test_framework.hpp"
#include "zeta_mobius.h"
#include "multivector.h"
#include <vector>
#include <array>
#include <random>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <span>
#include <numeric>

namespace moire_test {

// Test suite for ZetaMobius
class ZetaMobiusTestSuite : public TestSuite {
public:
    ZetaMobiusTestSuite() : TestSuite("ZetaMobius") {
        add_test("MultinomialCalculator Basic", [this]() { test_multinomial_calculator_basic(); });
        add_test("MultinomialCalculator Edge Cases", [this]() { test_multinomial_calculator_edge_cases(); });
        add_test("MobiusTransform Construction", [this]() { test_mobius_transform_construction(); });
        add_test("MobiusTransform Query Probabilities", [this]() { test_mobius_transform_query_probabilities(); });
        add_test("MobiusTransform Update Parameters", [this]() { test_mobius_transform_update_parameters(); });
        add_test("Support Vector Conversion", [this]() { test_support_vector_conversion(); });
        add_test("Mathematical Properties", [this]() { test_mathematical_properties(); });
        add_test("Edge Cases", [this]() { test_edge_cases(); });
        add_test("Performance", [this]() { test_performance(); });
        add_test("Create Multinomial Mobius Transform", [this]() { test_create_multinomial_mobius_transform(); });
        add_test("Probability Calculator Concept", [this]() { test_probability_calculator_concept(); });
        add_test("Large Scale Tests", [this]() { test_large_scale(); });
        add_test("Expected Data Tests", [this]() { test_expected_data(); });
        add_test("Mobius Transform Mathematical Tests", [this]() { test_mobius_transform_mathematical(); });
    }

private:
    void test_multinomial_calculator_basic() {
        std::vector<double> probs = {0.3, 0.4, 0.3};
        std::size_t max_n = 5;
        std::size_t u_mask = 0b111; // All categories included
        
        zeta_mobius::MultinomialCalculator calc;
        auto result = calc(probs, max_n, u_mask);
        
        // Debug logging
        std::cout << "DEBUG: probs = {0.3, 0.4, 0.3}, max_n = " << max_n << ", u_mask = " << u_mask << std::endl;
        std::cout << "DEBUG: result.size() = " << result.size() << std::endl;
        for (size_t i = 0; i < result.size(); ++i) {
            std::cout << "DEBUG: result[" << i << "] = " << result[i] << std::endl;
        }
        
        // Check that result has the correct size
        ASSERT_EQ(result.size(), max_n + 1);
        
        // Check that first element is 0.0 (not used in calculations)
        std::cout << "DEBUG: Expected result[0] = 0.0, actual = " << result[0] << std::endl;
        ASSERT_EQ(result[0], 0.0);
        
        // Check that second element is the probability that a single event is not in the support
        // With u_mask = 0b111 (all categories), result[1] should be 1.0 - (0.3 + 0.4 + 0.3) = 0.0
        std::cout << "DEBUG: Expected result[1] = 0.0, actual = " << result[1] << std::endl;
        ASSERT_EQ(result[1], 0.0);
        
        // Check that the result is non-empty and has the right size
        ASSERT_GT(result.size(), 0);
        
        // With the updated implementation, we can check that the calculation is more reasonable
        // The first element should be negative (as it's subtracted from)
        ASSERT_LE(result[0], 0.0);
    }

    void test_multinomial_calculator_edge_cases() {
        zeta_mobius::MultinomialCalculator calc;
        
        // Test with single category
        std::vector<double> single_prob = {1.0};
        auto result1 = calc(single_prob, 3, 0b1);
        ASSERT_EQ(result1.size(), 4);
        
        // Debug logging
        std::cout << "DEBUG: single_prob = {1.0}, u_mask = 0b1" << std::endl;
        for (size_t i = 0; i < result1.size(); ++i) {
            std::cout << "DEBUG: result1[" << i << "] = " << result1[i] << std::endl;
        }
        
        // With single category and u_mask = 0b1, result[0] should be 0.0 (not used)
        std::cout << "DEBUG: Expected result1[0] = 0.0, actual = " << result1[0] << std::endl;
        ASSERT_EQ(result1[0], 0.0);
        // Second element should be 1.0 - 1.0 = 0.0 (probability that single event is not in support)
        std::cout << "DEBUG: Expected result1[1] = 0.0, actual = " << result1[1] << std::endl;
        ASSERT_EQ(result1[1], 0.0);
        
        // Test with zero mask (no categories)
        std::vector<double> probs = {0.5, 0.5};
        auto result2 = calc(probs, 2, 0b00);
        ASSERT_EQ(result2.size(), 3);
        // With zero mask, no probabilities are subtracted, so result[0] is 0.0 (not used)
        ASSERT_EQ(result2[0], 0.0);
        // Second element should be 1.0 (no probabilities subtracted, so single event not in support = 1.0)
        ASSERT_EQ(result2[1], 1.0);
        
        // Test with partial mask
        auto result3 = calc(probs, 2, 0b01); // Only first category
        ASSERT_EQ(result3.size(), 3);
        // With partial mask (0b01), result[0] should be 0.0 (not used)
        ASSERT_EQ(result3[0], 0.0);
        // Second element should be 1.0 - 0.5 = 0.5 (probability that single event is not in support)
        ASSERT_EQ(result3[1], 0.5);
    }

    void test_mobius_transform_construction() {
        std::vector<double> probs = {0.3, 0.4, 0.3};
        std::size_t max_n = 5;
        
        auto transform = zeta_mobius::create_multinomial_mobius_transform(probs, max_n);
        
        // Test that we can query probabilities
        std::vector<int> support_vector = {1, 1, 1};
        auto prob_span = transform.query_probabilities(support_vector);
        
        ASSERT_GT(prob_span.size(), 0);
        
        // Check that probabilities are non-negative
        for (size_t i = 0; i < prob_span.size(); ++i) {
            ASSERT_GE(prob_span[i], 0.0);
        }
    }

    void test_mobius_transform_query_probabilities() {
        std::vector<double> probs = {0.5, 0.5};
        std::size_t max_n = 3;
        
        auto transform = zeta_mobius::create_multinomial_mobius_transform(probs, max_n);
        
        // Test different support vectors
        std::vector<std::vector<int>> support_vectors = {
            {1, 1},  // Both categories
            {1, 0},  // First category only
            {0, 1},  // Second category only
            {0, 0}   // No categories
        };
        
        for (const auto& support : support_vectors) {
            auto prob_span = transform.query_probabilities(support);
            
            // Check that we get the expected number of probabilities (inner dimension size is max_n + 1)
            ASSERT_EQ(prob_span.size(), max_n + 1);
            
            // Check that probabilities are non-negative
            for (size_t i = 0; i < prob_span.size(); ++i) {
                ASSERT_GE(prob_span[i], 0.0);
            }
        }
    }

    void test_mobius_transform_update_parameters() {
        std::vector<double> initial_probs = {0.3, 0.7};
        std::size_t initial_max_n = 3;
        
        auto transform = zeta_mobius::create_multinomial_mobius_transform(initial_probs, initial_max_n);
        
        // Query with initial parameters
        std::vector<int> support = {1, 1};
        auto initial_probs_span = transform.query_probabilities(support);
        
        // Update parameters
        std::vector<double> new_probs = {0.6, 0.4};
        std::size_t new_max_n = 4;
        transform.update_parameters(new_probs, new_max_n);
        
        // Query with new parameters
        auto new_probs_span = transform.query_probabilities(support);
        
        // Check that the new span has the correct size (inner dimension size is max_n + 1)
        ASSERT_EQ(new_probs_span.size(), new_max_n + 1);
        ASSERT_GT(new_probs_span.size(), initial_probs_span.size());
    }

    void test_support_vector_conversion() {
        // This tests the internal conversion functions indirectly
        std::vector<double> probs = {0.25, 0.25, 0.25, 0.25};
        std::size_t max_n = 2;
        
        auto transform = zeta_mobius::create_multinomial_mobius_transform(probs, max_n);
        
        // Test various support vector patterns
        std::vector<std::vector<int>> test_vectors = {
            {1, 0, 0, 0},
            {0, 1, 0, 0},
            {0, 0, 1, 0},
            {0, 0, 0, 1},
            {1, 1, 0, 0},
            {1, 0, 1, 0},
            {0, 1, 1, 0},
            {1, 1, 1, 1}
        };
        
        for (const auto& support : test_vectors) {
            auto prob_span = transform.query_probabilities(support);
            
            // Each query should return valid probabilities (inner dimension size is max_n + 1)
            ASSERT_EQ(prob_span.size(), max_n + 1);
            for (size_t i = 0; i < prob_span.size(); ++i) {
                ASSERT_TRUE(std::isfinite(prob_span[i]));
            }
        }
    }

    void test_mathematical_properties() {
        std::vector<double> probs = {0.4, 0.3, 0.3};
        std::size_t max_n = 4;
        
        auto transform = zeta_mobius::create_multinomial_mobius_transform(probs, max_n);
        
        // Test different support vectors
        std::vector<std::vector<int>> support_vectors = {
            {1, 0, 0},
            {0, 1, 0},
            {0, 0, 1},
            {1, 1, 0},
            {1, 0, 1},
            {0, 1, 1},
            {1, 1, 1}
        };
        
        for (const auto& support : support_vectors) {
            auto prob_span = transform.query_probabilities(support);
            
            // Check that we get the expected size
            ASSERT_EQ(prob_span.size(), max_n + 1);
            
            // Check that probabilities are non-negative
            for (size_t i = 0; i < prob_span.size(); ++i) {
                ASSERT_GE(prob_span[i], 0.0);
            }
            
            // The first element (index 0) represents "no events" which should be 0
            // for non-empty support vectors
            if (std::any_of(support.begin(), support.end(), [](int x) { return x != 0; })) {
                ASSERT_EQ(prob_span[0], 0.0);
            }
        }
    }

    void test_edge_cases() {
        // Test with very small probabilities
        std::vector<double> small_probs = {0.001, 0.999};
        std::size_t max_n = 2;
        
        auto transform1 = zeta_mobius::create_multinomial_mobius_transform(small_probs, max_n);
        std::vector<int> support = {1, 1};
        auto prob_span1 = transform1.query_probabilities(support);
        
        ASSERT_EQ(prob_span1.size(), max_n + 1);
        for (size_t i = 0; i < prob_span1.size(); ++i) {
            ASSERT_GE(prob_span1[i], 0.0);
        }
        
        // Test with equal probabilities
        std::vector<double> equal_probs = {0.5, 0.5};
        auto transform2 = zeta_mobius::create_multinomial_mobius_transform(equal_probs, max_n);
        auto prob_span2 = transform2.query_probabilities(support);
        
        ASSERT_EQ(prob_span2.size(), max_n + 1);
        
        // Test with single category
        std::vector<double> single_prob = {1.0};
        auto transform3 = zeta_mobius::create_multinomial_mobius_transform(single_prob, 1);
        std::vector<int> single_support = {1};
        auto prob_span3 = transform3.query_probabilities(single_support);
        
        ASSERT_EQ(prob_span3.size(), 2);
        // First element should be 0 (no events for non-empty support)
        ASSERT_EQ(prob_span3[0], 0.0);
        // Second element should be non-negative
        ASSERT_GE(prob_span3[1], 0.0);
    }

    void test_performance() {
        // Test with larger parameters
        std::vector<double> probs(8, 0.125); // 8 categories with equal probability
        std::size_t max_n = 10;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        auto transform = zeta_mobius::create_multinomial_mobius_transform(probs, max_n);
        
        // Query multiple support vectors
        std::vector<std::vector<int>> support_vectors = {
            {1, 0, 0, 0, 0, 0, 0, 0},
            {1, 1, 0, 0, 0, 0, 0, 0},
            {1, 1, 1, 0, 0, 0, 0, 0},
            {1, 1, 1, 1, 0, 0, 0, 0},
            {1, 1, 1, 1, 1, 0, 0, 0},
            {1, 1, 1, 1, 1, 1, 0, 0},
            {1, 1, 1, 1, 1, 1, 1, 0},
            {1, 1, 1, 1, 1, 1, 1, 1}
        };
        
        for (const auto& support : support_vectors) {
            auto prob_span = transform.query_probabilities(support);
            ASSERT_EQ(prob_span.size(), max_n + 1);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        
        // Performance should be reasonable (less than 1 second for this test)
        ASSERT_LT(duration, 1000);
    }

    void test_create_multinomial_mobius_transform() {
        std::vector<double> probs = {0.2, 0.3, 0.5};
        std::size_t max_n = 3;
        
        auto transform = zeta_mobius::create_multinomial_mobius_transform(probs, max_n);
        
        // Test that the transform works correctly
        std::vector<int> support = {1, 1, 1};
        auto prob_span = transform.query_probabilities(support);
        
        ASSERT_EQ(prob_span.size(), max_n + 1);
        
        // Test update functionality
        std::vector<double> new_probs = {0.1, 0.4, 0.5};
        transform.update_parameters(new_probs, max_n);
        
        auto new_prob_span = transform.query_probabilities(support);
        
        // Debug logging
        std::cout << "DEBUG: new_max_n = " << max_n << ", new_prob_span.size() = " << new_prob_span.size() << std::endl;
        for (size_t i = 0; i < new_prob_span.size(); ++i) {
            std::cout << "DEBUG: new_prob_span[" << i << "] = " << new_prob_span[i] << std::endl;
        }
        
        ASSERT_EQ(new_prob_span.size(), max_n + 1);
    }

    void test_probability_calculator_concept() {
        // Test that MultinomialCalculator satisfies the ProbabilityCalculator concept
        zeta_mobius::MultinomialCalculator calc;
        std::vector<double> probs = {0.5, 0.5};
        std::size_t max_n = 2;
        std::size_t u_mask = 0b11;
        
        // This should compile and work if the concept is satisfied
        auto result = calc(probs, max_n, u_mask);
        ASSERT_EQ(result.size(), max_n + 1);
    }

    void test_large_scale() {
        // Test with larger scale parameters
        const size_t num_categories = 6;
        const size_t max_n = 8;
        
        std::vector<double> probs(num_categories);
        std::iota(probs.begin(), probs.end(), 1.0); // 1, 2, 3, 4, 5, 6
        double sum = std::accumulate(probs.begin(), probs.end(), 0.0);
        for (auto& p : probs) {
            p /= sum; // Normalize to sum to 1
        }
        
        auto transform = zeta_mobius::create_multinomial_mobius_transform(probs, max_n);
        
        // Test various support patterns
        std::vector<std::vector<int>> patterns = {
            {1, 0, 0, 0, 0, 0},
            {0, 1, 0, 0, 0, 0},
            {1, 1, 0, 0, 0, 0},
            {1, 0, 1, 0, 0, 0},
            {1, 1, 1, 0, 0, 0},
            {1, 1, 1, 1, 0, 0},
            {1, 1, 1, 1, 1, 0},
            {1, 1, 1, 1, 1, 1}
        };
        
        for (const auto& pattern : patterns) {
            auto prob_span = transform.query_probabilities(pattern);
            
            ASSERT_EQ(prob_span.size(), max_n + 1);
            
            // Check that probabilities are non-negative
            for (size_t i = 0; i < prob_span.size(); ++i) {
                ASSERT_GE(prob_span[i], 0.0);
            }
        }
        
        // Test parameter updates
        std::vector<double> new_probs(num_categories, 1.0 / num_categories);
        transform.update_parameters(new_probs, max_n + 1);
        
        auto new_prob_span = transform.query_probabilities(std::span<const int>({1, 1, 1, 1, 1, 1}));
        
        // Debug logging
        std::cout << "DEBUG: max_n = " << max_n << ", new_prob_span.size() = " << new_prob_span.size() << std::endl;
        for (size_t i = 0; i < new_prob_span.size(); ++i) {
            std::cout << "DEBUG: new_prob_span[" << i << "] = " << new_prob_span[i] << std::endl;
        }
        
        ASSERT_EQ(new_prob_span.size(), max_n + 2);
    }

    void test_expected_data() {
        zeta_mobius::MultinomialCalculator calc;
        
        // Test 1: Single category with probability 0.5, not in support
        // Expected: P(not in support) = 1.0 - 0.5 = 0.5
        // For n events: P(not in support)^n = 0.5^n
        std::vector<double> probs1 = {0.5};
        auto result1 = calc(probs1, 3, 0b0); // No categories in support
        ASSERT_EQ(result1.size(), 4);
        ASSERT_EQ(result1[0], 0.0); // Not used
        ASSERT_EQ(result1[1], 1.0); // 1.0 - 0.0 = 1.0 (no categories in support)
        ASSERT_EQ(result1[2], 1.0); // 1.0^2 = 1.0
        ASSERT_EQ(result1[3], 1.0); // 1.0^3 = 1.0
        
        // Test 2: Single category with probability 0.5, in support
        // Expected: P(not in support) = 1.0 - 0.5 = 0.5
        // For n events: P(not in support)^n = 0.5^n
        auto result2 = calc(probs1, 3, 0b1); // Category in support
        ASSERT_EQ(result2.size(), 4);
        ASSERT_EQ(result2[0], 0.0); // Not used
        ASSERT_TRUE(std::abs(result2[1] - 0.5) < 1e-10); // 1.0 - 0.5 = 0.5
        ASSERT_TRUE(std::abs(result2[2] - 0.25) < 1e-10); // 0.5^2 = 0.25
        ASSERT_TRUE(std::abs(result2[3] - 0.125) < 1e-10); // 0.5^3 = 0.125
        
        // Test 3: Two categories with equal probabilities, one in support
        // Categories: [0.4, 0.6], support: [1, 0] (first category in support)
        // Expected: P(not in support) = 1.0 - 0.4 = 0.6
        std::vector<double> probs3 = {0.4, 0.6};
        auto result3 = calc(probs3, 2, 0b01); // First category in support
        ASSERT_EQ(result3.size(), 3);
        ASSERT_EQ(result3[0], 0.0); // Not used
        ASSERT_TRUE(std::abs(result3[1] - 0.6) < 1e-10); // 1.0 - 0.4 = 0.6
        ASSERT_TRUE(std::abs(result3[2] - 0.36) < 1e-10); // 0.6^2 = 0.36
        
        // Test 4: Two categories with equal probabilities, both in support
        // Categories: [0.4, 0.6], support: [1, 1] (both categories in support)
        // Expected: P(not in support) = 1.0 - 0.4 - 0.6 = 0.0
        auto result4 = calc(probs3, 2, 0b11); // Both categories in support
        ASSERT_EQ(result4.size(), 3);
        ASSERT_EQ(result4[0], 0.0); // Not used
        ASSERT_EQ(result4[1], 0.0); // 1.0 - 0.4 - 0.6 = 0.0
        ASSERT_EQ(result4[2], 0.0); // 0.0^2 = 0.0
        
        // Test 5: Three categories with different probabilities, partial support
        // Categories: [0.2, 0.3, 0.5], support: [1, 0, 1] (first and third in support)
        // Expected: P(not in support) = 1.0 - 0.2 - 0.5 = 0.3
        std::vector<double> probs5 = {0.2, 0.3, 0.5};
        auto result5 = calc(probs5, 3, 0b101); // First and third categories in support
        ASSERT_EQ(result5.size(), 4);
        ASSERT_EQ(result5[0], 0.0); // Not used
        
        // Debug output to see actual values
        std::cout << "DEBUG: result5[1] = " << result5[1] << ", expected = 0.3" << std::endl;
        std::cout << "DEBUG: result5[2] = " << result5[2] << ", expected = 0.09" << std::endl;
        std::cout << "DEBUG: result5[3] = " << result5[3] << ", expected = 0.027" << std::endl;
        
        // Use approximate equality for floating-point comparisons
        ASSERT_TRUE(std::abs(result5[1] - 0.3) < 1e-10); // 1.0 - 0.2 - 0.5 = 0.3
        ASSERT_TRUE(std::abs(result5[2] - 0.09) < 1e-10); // 0.3^2 = 0.09
        ASSERT_TRUE(std::abs(result5[3] - 0.027) < 1e-10); // 0.3^3 = 0.027
        
        // Test 6: Edge case with very small probabilities
        // Categories: [0.1, 0.9], support: [1, 0] (first category in support)
        // Expected: P(not in support) = 1.0 - 0.1 = 0.9
        std::vector<double> probs6 = {0.1, 0.9};
        auto result6 = calc(probs6, 2, 0b01); // First category in support
        ASSERT_EQ(result6.size(), 3);
        ASSERT_EQ(result6[0], 0.0); // Not used
        ASSERT_TRUE(std::abs(result6[1] - 0.9) < 1e-10); // 1.0 - 0.1 = 0.9
        ASSERT_TRUE(std::abs(result6[2] - 0.81) < 1e-10); // 0.9^2 = 0.81
        
        // Test 7: Edge case with probabilities that sum to 1.0
        // Categories: [0.3, 0.7], support: [1, 1] (both categories in support)
        // Expected: P(not in support) = 1.0 - 0.3 - 0.7 = 0.0
        std::vector<double> probs7 = {0.3, 0.7};
        auto result7 = calc(probs7, 2, 0b11); // Both categories in support
        ASSERT_EQ(result7.size(), 3);
        ASSERT_EQ(result7[0], 0.0); // Not used
        ASSERT_EQ(result7[1], 0.0); // 1.0 - 0.3 - 0.7 = 0.0
        ASSERT_EQ(result7[2], 0.0); // 0.0^2 = 0.0
    }

    void test_mobius_transform_mathematical() {
        // Test 1: Single category support vector - inclusion-exclusion principle
        // The Mobius transform implements inclusion-exclusion for probability calculations
        std::vector<double> probs = {0.3, 0.4, 0.3};
        std::size_t max_n = 3;
        auto transform = zeta_mobius::create_multinomial_mobius_transform(probs, max_n);
        
        std::vector<int> single_support = {1, 0, 0}; // Only first category in support
        auto single_span = transform.query_probabilities(single_support);
        ASSERT_EQ(single_span.size(), max_n + 1);
        
        // Debug output to see actual values
        std::cout << "DEBUG: Single support test - actual values:" << std::endl;
        for (size_t i = 0; i < single_span.size(); ++i) {
            std::cout << "DEBUG: single_span[" << i << "] = " << single_span[i] << std::endl;
        }
        
        // The Mobius transform calculates probabilities using inclusion-exclusion
        // For a single category support, this calculates the probability of events IN the support
        // Expected: P(in support) = 0.3 (the probability of the first category)
        ASSERT_EQ(single_span[0], 0.0); // Not used
        ASSERT_TRUE(std::abs(single_span[1] - 0.3) < 1e-10); // P(in support) = 0.3
        ASSERT_TRUE(std::abs(single_span[2] - 0.09) < 1e-10); // 0.3^2 = 0.09
        ASSERT_TRUE(std::abs(single_span[3] - 0.027) < 1e-10); // 0.3^3 = 0.027
        
        // Test 2: Full support vector - inclusion-exclusion with all categories
        // When all categories are in support, the Mobius transform should handle
        // the inclusion-exclusion principle correctly
        std::vector<int> full_support = {1, 1, 1};
        auto full_span = transform.query_probabilities(full_support);
        ASSERT_EQ(full_span.size(), max_n + 1);
        
        // With all categories in support, the inclusion-exclusion principle
        // should result in probabilities of 1.0 (all events are in the support)
        for (size_t i = 0; i < full_span.size(); ++i) {
            if (i == 0) {
                ASSERT_EQ(full_span[i], 0.0); // Not used
            } else {
                ASSERT_TRUE(std::abs(full_span[i] - 1.0) < 1e-10); // 1.0^i = 1.0
            }
        }
        
        // Test 3: Partial support vector - inclusion-exclusion with multiple categories
        // The Mobius transform handles inclusion-exclusion for overlapping events
        // Categories: [0.2, 0.3, 0.5], support: [1, 0, 1] (first and third in support)
        std::vector<double> probs3 = {0.2, 0.3, 0.5};
        auto transform3 = zeta_mobius::create_multinomial_mobius_transform(probs3, 3);
        std::vector<int> partial_support = {1, 0, 1};
        auto partial_span = transform3.query_probabilities(partial_support);
        ASSERT_EQ(partial_span.size(), 4);
        
        // The Mobius transform applies inclusion-exclusion principle:
        // P(in support) = P(first) + P(third) = 0.2 + 0.5 = 0.7
        // This avoids double-counting overlapping events
        ASSERT_EQ(partial_span[0], 0.0); // Not used
        ASSERT_TRUE(std::abs(partial_span[1] - 0.7) < 1e-10); // 0.2 + 0.5 = 0.7
        ASSERT_TRUE(std::abs(partial_span[2] - 0.49) < 1e-10); // 0.7^2 = 0.49
        ASSERT_TRUE(std::abs(partial_span[3] - 0.343) < 1e-10); // 0.7^3 = 0.343
        
        // Test 4: Single category support
        std::vector<double> probs4 = {0.6, 0.4};
        auto transform4 = zeta_mobius::create_multinomial_mobius_transform(probs4, 2);
        std::vector<int> single_support_2 = {1, 0};
        auto single_span_2 = transform4.query_probabilities(single_support_2);
        ASSERT_EQ(single_span_2.size(), 3);
        
        // Expected: P(in support) = 0.6 (the probability of the first category)
        ASSERT_EQ(single_span_2[0], 0.0); // Not used
        ASSERT_TRUE(std::abs(single_span_2[1] - 0.6) < 1e-10); // P(in support) = 0.6
        ASSERT_TRUE(std::abs(single_span_2[2] - 0.36) < 1e-10); // 0.6^2 = 0.36
        
        // Test 5: Mobius transform inclusion-exclusion properties
        // The Mobius transform should satisfy inclusion-exclusion mathematical properties
        std::vector<double> probs5 = {0.25, 0.25, 0.25, 0.25};
        auto transform5 = zeta_mobius::create_multinomial_mobius_transform(probs5, 2);
        
        // Test different support vectors to verify inclusion-exclusion principle
        std::vector<std::vector<int>> test_supports = {
            {1, 0, 0, 0}, // First category only
            {0, 1, 0, 0}, // Second category only
            {0, 0, 1, 0}, // Third category only
            {1, 1, 0, 0}, // First two categories (tests inclusion-exclusion)
            {1, 1, 1, 1}  // All categories (tests full inclusion-exclusion)
        };
        
        for (const auto& support : test_supports) {
            auto span = transform5.query_probabilities(support);
            ASSERT_EQ(span.size(), 3);
            
            // All probabilities should be non-negative (inclusion-exclusion property)
            for (size_t i = 0; i < span.size(); ++i) {
                ASSERT_GE(span[i], 0.0);
            }
            
            // The first element should be 0.0 (not used)
            ASSERT_EQ(span[0], 0.0);
            
            // The probabilities should be monotonically decreasing (inclusion-exclusion property)
            // This reflects the mathematical structure of the Mobius transform
            if (span[1] > 0.0) {
                ASSERT_LE(span[2], span[1]);
            }
        }
        
        // Test 6: Parameter update consistency with inclusion-exclusion
        // The Mobius transform should maintain inclusion-exclusion properties after updates
        std::vector<double> initial_probs = {0.5, 0.5};
        auto transform6 = zeta_mobius::create_multinomial_mobius_transform(initial_probs, 2);
        std::vector<int> support = {1, 0};
        auto initial_span = transform6.query_probabilities(support);
        
        // Update parameters - the inclusion-exclusion principle should still apply
        std::vector<double> new_probs = {0.3, 0.7};
        transform6.update_parameters(new_probs, 2);
        auto updated_span = transform6.query_probabilities(support);
        
        // Both should have the same size (inclusion-exclusion structure preserved)
        ASSERT_EQ(initial_span.size(), updated_span.size());
        
        // The updated probabilities should reflect the new parameters using inclusion-exclusion
        // With new_probs = {0.3, 0.7} and support = {1, 0}
        // Expected: P(in support) = 0.3 (inclusion-exclusion principle)
        ASSERT_TRUE(std::abs(updated_span[1] - 0.3) < 1e-10); // P(in support) = 0.3
        ASSERT_TRUE(std::abs(updated_span[2] - 0.09) < 1e-10); // 0.3^2 = 0.09
        
        // Test 7: Edge case with very small probabilities - inclusion-exclusion robustness
        // The Mobius transform should handle extreme probability values correctly
        std::vector<double> small_probs = {0.01, 0.99};
        auto transform7 = zeta_mobius::create_multinomial_mobius_transform(small_probs, 2);
        std::vector<int> small_support = {1, 0};
        auto small_span = transform7.query_probabilities(small_support);
        
        // Expected: P(in support) = 0.01 (inclusion-exclusion principle)
        ASSERT_TRUE(std::abs(small_span[1] - 0.01) < 1e-10); // P(in support) = 0.01
        ASSERT_TRUE(std::abs(small_span[2] - 0.0001) < 1e-10); // 0.01^2 = 0.0001
        
        // Test 8: Mathematical consistency across different max_n values
        // The inclusion-exclusion principle should be consistent regardless of max_n
        std::vector<double> probs8 = {0.4, 0.6};
        auto transform8_small = zeta_mobius::create_multinomial_mobius_transform(probs8, 2);
        auto transform8_large = zeta_mobius::create_multinomial_mobius_transform(probs8, 4);
        
        std::vector<int> support8 = {1, 0};
        auto span_small = transform8_small.query_probabilities(support8);
        auto span_large = transform8_large.query_probabilities(support8);
        
        // The first few elements should be identical (inclusion-exclusion consistency)
        for (size_t i = 0; i < span_small.size(); ++i) {
            ASSERT_TRUE(std::abs(span_small[i] - span_large[i]) < 1e-10);
        }
        
        // The larger span should have additional elements (inclusion-exclusion structure)
        ASSERT_GT(span_large.size(), span_small.size());
    }
};

} // namespace moire_test

int main() {
    auto suite = std::make_unique<moire_test::ZetaMobiusTestSuite>();
    moire_test::TestRegistry::register_suite(std::move(suite));
    moire_test::TestRegistry::run_all_tests();
    return 0;
}
