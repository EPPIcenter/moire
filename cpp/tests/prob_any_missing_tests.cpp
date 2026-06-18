#include "prob_any_missing.h"
#include "prob_any_missing_cache.h"
#include "../prob_any_missing_alternatives.h"
#include "test_framework.hpp"

#include <cmath>
#include <vector>

namespace moire_test {

class ProbAnyMissingTestSuite : public TestSuite {
public:
    ProbAnyMissingTestSuite() : TestSuite("ProbAnyMissing") {
        add_test("Vectorized matches combination reference", [this]() {
            test_vectorized_matches_combination();
        });
        add_test("Operator matches vectorized at coi index", [this]() {
            test_operator_matches_vectorized();
        });
        add_test("Impossible coi returns certainty", [this]() {
            test_impossible_coi_returns_one();
        });
        add_test("PamVectorCache returns stored vectorized result", [this]() {
            test_pam_vector_cache();
        });
        add_test("Quantized simplex buckets nearby q", [this]() {
            test_quantized_simplex_buckets();
        });
    }

private:
    static constexpr double kTol = 1e-5;

    static void assert_near(double actual, double expected, const char* label)
    {
        if (std::fabs(actual - expected) > kTol) {
            throw std::runtime_error(
                std::string(label) + ": expected " + std::to_string(expected) + ", got "
                + std::to_string(actual));
        }
    }

    void test_vectorized_matches_combination()
    {
        probAnyMissingFunctor production;
        prob_any_missing_bench::ReferenceFunctor reference;

        const std::vector<std::vector<float>> fixtures = {
            {0.5f, 0.5f},
            {0.2f, 0.3f, 0.5f},
            {0.1f, 0.2f, 0.3f, 0.4f},
        };
        const std::vector<unsigned> coi_values = {1u, 2u, 5u};

        for (const auto& event_probs : fixtures) {
            for (unsigned coi : coi_values) {
                if (coi < event_probs.size()) continue;
                const auto gray = production.vectorized(event_probs, coi);
                const auto comb = reference.vectorized_combination(event_probs, coi);
                if (gray.size() != comb.size()) {
                    throw std::runtime_error("size mismatch for coi=" + std::to_string(coi));
                }
                for (std::size_t i = 0; i < gray.size(); ++i) {
                    assert_near(gray[i], comb[i], "gray vs combination");
                }
            }
        }
    }

    void test_operator_matches_vectorized()
    {
        probAnyMissingFunctor functor;
        const std::vector<float> event_probs = {0.25f, 0.25f, 0.5f};

        for (int coi = 1; coi <= 6; ++coi) {
            const double scalar = functor(event_probs, coi);
            const auto vec = functor.vectorized(event_probs, static_cast<unsigned>(coi));
            assert_near(scalar, vec[static_cast<std::size_t>(coi - 1)], "operator vs vectorized");
        }
    }

    void test_impossible_coi_returns_one()
    {
        probAnyMissingFunctor functor;
        const std::vector<float> event_probs = {0.4f, 0.6f};
        const auto result = functor.vectorized(event_probs, 1u);
        ASSERT_EQ(result.size(), 1u);
        assert_near(result[0], 1.0, "impossible coi");
    }

    void test_pam_vector_cache()
    {
        probAnyMissingFunctor functor;
        PamVectorCache cache;
        const std::vector<float> event_probs = {0.25f, 0.25f, 0.5f};
        const unsigned coi = 4u;
        const auto expected = functor.vectorized(event_probs, 1u, coi);

        ASSERT_TRUE(cache.lookup(event_probs, 1u, coi) == nullptr);
        const PamCachedVectors& stored =
            cache.store_and_get(event_probs, 1u, coi, expected);
        ASSERT_EQ(stored.pam.size(), expected.size());
        for (std::size_t i = 0; i < expected.size(); ++i) {
            assert_near(stored.pam[i], expected[i], "cache round-trip");
            assert_near(static_cast<double>(stored.log_one_minus_pam[i]),
                        std::log1p(-expected[i]), "log(1-pam)");
        }
        const PamCachedVectors* hit = cache.lookup(event_probs, 1u, coi);
        ASSERT_TRUE(hit != nullptr);
        ASSERT_EQ(hit->pam.size(), expected.size());
    }

    void test_quantized_simplex_buckets()
    {
        const std::vector<float> q0 = {0.25001f, 0.24999f, 0.50000f};
        const std::vector<float> q1 = {0.25004f, 0.24996f, 0.50000f};
        const std::vector<float> key0 = pam_cache::quantize_simplex(q0, 1e-4f);
        const std::vector<float> key1 = pam_cache::quantize_simplex(q1, 1e-4f);
        ASSERT_EQ(key0.size(), key1.size());
        for (std::size_t i = 0; i < key0.size(); ++i) {
            assert_near(key0[i], key1[i], "quantized bucket");
        }

        PamVectorCache cache;
        probAnyMissingFunctor functor;
        const unsigned coi = 4u;
        const auto exact0 = functor.vectorized(q0, 1u, coi);
        cache.store_and_get(q0, 1u, coi, exact0);
        ASSERT_TRUE(cache.lookup(q1, 1u, coi) != nullptr);
    }
};

} // namespace moire_test

int main()
{
    moire_test::ProbAnyMissingTestSuite suite;
    suite.run_tests();
    return suite.all_passed() ? 0 : 1;
}
