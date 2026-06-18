#include "test_framework.hpp"
#include "binary_observation_model.h"
#include "count_poisson_observation_model.h"

#include <cmath>
#include <vector>

namespace moire_test {

class ObservationModelTestSuite : public TestSuite {
public:
    ObservationModelTestSuite() : TestSuite("ObservationModel") {
        add_test("Binary perfect match two alleles", [this]() {
            test_binary_perfect_match_two_alleles();
        });
        add_test("Binary all true negative", [this]() { test_binary_all_true_negative(); });
        add_test("Binary false positive", [this]() { test_binary_false_positive(); });
        add_test("Binary false negative", [this]() { test_binary_false_negative(); });
        add_test("Binary mixed tp tn fp fn", [this]() { test_binary_mixed_tp_tn_fp_fn(); });
        add_test("Count absent noise only", [this]() { test_count_absent_noise_only(); });
        add_test("Count present signal", [this]() { test_count_present_signal(); });
        add_test("Count dropout mixture", [this]() { test_count_dropout_mixture(); });
        add_test("Binary latent padding ignored", [this]() {
            test_binary_latent_padding_ignored();
        });
        add_test("Count latent padding ignored", [this]() {
            test_count_latent_padding_ignored();
        });
    }

private:
    static constexpr float kMaxEpsNeg = 2.0f;
    static constexpr float kMaxEpsPos = 2.0f;
    static constexpr float kFixtureTolerance = 1e-5f;

    void assert_log_likelihood_near(float actual, float expected)
    {
        if (std::abs(actual - expected) > kFixtureTolerance) {
            throw std::runtime_error(
                "log-likelihood mismatch: expected " + std::to_string(expected) + ", got "
                + std::to_string(actual));
        }
    }

    void test_binary_perfect_match_two_alleles()
    {
        BinaryObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent = {0, 1, -1};
        std::vector<int> observed = {1, 1};

        const float actual =
            model.log_likelihood(latent, observed, 0.1f, 0.2f);
        assert_log_likelihood_near(actual, -0.210721f);
    }

    void test_binary_all_true_negative()
    {
        BinaryObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent = {-1};
        std::vector<int> observed = {0, 0};

        const float actual =
            model.log_likelihood(latent, observed, 0.05f, 0.15f);
        assert_log_likelihood_near(actual, -0.3250379f);
    }

    void test_binary_false_positive()
    {
        BinaryObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent = {-1};
        std::vector<int> observed = {1, 0};

        const float actual =
            model.log_likelihood(latent, observed, 0.1f, 0.1f);
        assert_log_likelihood_near(actual, -2.407946f);
    }

    void test_binary_false_negative()
    {
        BinaryObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent = {0, -1};
        std::vector<int> observed = {0, 0};

        const float actual =
            model.log_likelihood(latent, observed, 0.2f, 0.3f);
        assert_log_likelihood_near(actual, -1.966113f);
    }

    void test_binary_mixed_tp_tn_fp_fn()
    {
        BinaryObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent = {0, 2, -1};
        std::vector<int> observed = {1, 0, 1};

        const float actual =
            model.log_likelihood(latent, observed, 0.1f, 0.1f);
        assert_log_likelihood_near(actual, -0.2069786f);
    }

    void test_count_absent_noise_only()
    {
        CountPoissonObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent = {-1};
        std::vector<int> observed = {0, 2, 1};

        const float actual =
            model.log_likelihood(latent, observed, 0.1f, 0.2f);
        assert_log_likelihood_near(actual, -7.137856f);
    }

    void test_count_present_signal()
    {
        CountPoissonObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent = {0, 2, -1};
        std::vector<int> observed = {5, 1, 10};

        const float actual =
            model.log_likelihood(latent, observed, 0.05f, 0.1f);
        assert_log_likelihood_near(actual, -6.300591f);
    }

    void test_count_dropout_mixture()
    {
        CountPoissonObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent = {1, -1};
        std::vector<int> observed = {0, 8};

        const float actual =
            model.log_likelihood(latent, observed, 0.3f, 0.15f);
        assert_log_likelihood_near(actual, -3.485359f);
    }

    void test_binary_latent_padding_ignored()
    {
        BinaryObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent_short = {0, 1, -1};
        std::vector<int> latent_padded = {0, 1, -1, -1, -1};
        std::vector<int> observed = {1, 1};

        const float short_ll =
            model.log_likelihood(latent_short, observed, 0.1f, 0.2f);
        const float padded_ll =
            model.log_likelihood(latent_padded, observed, 0.1f, 0.2f);
        assert_log_likelihood_near(short_ll, -0.210721f);
        assert_log_likelihood_near(padded_ll, short_ll);
    }

    void test_count_latent_padding_ignored()
    {
        CountPoissonObservationModel model(kMaxEpsNeg, kMaxEpsPos);
        std::vector<int> latent_short = {0, 2, -1};
        std::vector<int> latent_padded = {0, 2, -1, -1};
        std::vector<int> observed = {5, 1, 10};

        const float short_ll =
            model.log_likelihood(latent_short, observed, 0.05f, 0.1f);
        const float padded_ll =
            model.log_likelihood(latent_padded, observed, 0.05f, 0.1f);
        assert_log_likelihood_near(short_ll, -6.300591f);
        assert_log_likelihood_near(padded_ll, short_ll);
    }
};

}  // namespace moire_test

int main()
{
    moire_test::ObservationModelTestSuite suite;
    suite.run_tests();
    return suite.all_passed() ? 0 : 1;
}
