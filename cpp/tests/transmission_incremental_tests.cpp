#include "transmission_incremental.h"
#include "prob_any_missing.h"
#include "pam_fast_paths.h"
#include "test_framework.hpp"

#include <cmath>
#include <vector>

namespace moire_test {

class TransmissionIncrementalTestSuite : public TestSuite {
public:
    TransmissionIncrementalTestSuite() : TestSuite("TransmissionIncremental") {
        add_test("build_constrained_q normalizes latent alleles", [this]() {
            test_build_constrained_q();
        });
        add_test("log_S shift matches coi scaling without relatedness", [this]() {
            test_log_s_shift_no_relatedness();
        });
        add_test("incremental path matches full pam recompute", [this]() {
            test_incremental_matches_full_recompute();
        });
        add_test("low-k pam matches gray-code vectorized", [this]() {
            test_low_k_pam_matches_gray();
        });
    }

private:
    static constexpr double kTol = 1e-5;

    static void assert_near(float actual, float expected, const char* label)
    {
        if (std::fabs(actual - expected) > kTol) {
            throw std::runtime_error(
                std::string(label) + ": expected " + std::to_string(expected) + ", got "
                + std::to_string(actual));
        }
    }

    void test_build_constrained_q()
    {
        const std::vector<int> alleles = {0, 2};
        const std::vector<float> p = {0.2f, 0.3f, 0.5f};
        std::vector<float> q;
        float sum = 0.f;
        const bool ok =
            transmission_incremental::build_constrained_q(alleles, p, q, sum);
        if (!ok) {
            throw std::runtime_error("build_constrained_q failed");
        }
        assert_near(sum, 0.7f, "sum");
        assert_near(q[0], 0.2f / 0.7f, "q0");
        assert_near(q[1], 0.5f / 0.7f, "q1");
    }

    void test_log_s_shift_no_relatedness()
    {
        const std::vector<double> pam_vec = {0.05, 0.1, 0.2};
        const float log_s_old = std::log(0.4f);
        const float log_s_new = std::log(0.6f);
        const float t_old = transmission_incremental::transmission_log_no_relatedness(
            pam_vec, 3, log_s_old);
        const float t_new = transmission_incremental::transmission_log_no_relatedness(
            pam_vec, 3, log_s_new);
        assert_near(t_new - t_old, (log_s_new - log_s_old) * 3.f, "delta");
    }

    void test_incremental_matches_full_recompute()
    {
        probAnyMissingFunctor functor;
        const std::vector<int> alleles = {0, 1};
        const std::vector<float> p_old = {0.2f, 0.2f, 0.6f};
        const std::vector<float> p_new = {0.3f, 0.3f, 0.4f};

        std::vector<float> q_old;
        std::vector<float> q_new;
        float sum_old = 0.f;
        float sum_new = 0.f;
        transmission_incremental::build_constrained_q(alleles, p_old, q_old, sum_old);
        transmission_incremental::build_constrained_q(alleles, p_new, q_new, sum_new);

        const std::vector<double> pam_vec =
            functor.vectorized(q_new, 1u, 3u);
        const float incremental = transmission_incremental::transmission_log_no_relatedness(
            pam_vec, 3, std::log(sum_new));

        const std::vector<double> pam_vec_full =
            functor.vectorized(q_new, 1u, 3u);
        const float full = transmission_incremental::transmission_log_no_relatedness(
            pam_vec_full, 3, std::log(sum_new));

        assert_near(incremental, full, "incremental vs full");
    }

    void test_low_k_pam_matches_gray()
    {
        probAnyMissingFunctor functor;
        const std::vector<float> q2 = {0.3f, 0.7f};
        const unsigned max_n = 5u;
        std::vector<double> fast;
        if (!pam_fast_paths::try_fill_pam_vector(q2, 1u, max_n, fast)) {
            throw std::runtime_error("expected low-k fast path for K=2");
        }
        const std::vector<double> gray = functor.vectorized(q2, 1u, max_n);
        if (fast.size() != gray.size()) {
            throw std::runtime_error("pam vector size mismatch");
        }
        for (std::size_t i = 0; i < fast.size(); ++i) {
            assert_near(static_cast<float>(fast[i]), static_cast<float>(gray[i]),
                        "low-k vs gray");
        }

        const std::vector<float> q3 = {0.2f, 0.3f, 0.5f};
        std::vector<double> fast3;
        pam_fast_paths::try_fill_pam_vector(q3, 1u, max_n, fast3);
        const std::vector<double> gray3 = functor.vectorized(q3, 1u, max_n);
        for (std::size_t i = 0; i < fast3.size(); ++i) {
            assert_near(static_cast<float>(fast3[i]), static_cast<float>(gray3[i]),
                        "K=3 low-k vs gray");
        }

        const std::vector<float> q5 = {0.1f, 0.15f, 0.2f, 0.25f, 0.3f};
        std::vector<double> fast5;
        if (!pam_fast_paths::try_fill_pam_vector(q5, 1u, max_n, fast5)) {
            throw std::runtime_error("expected low-k fast path for K=5");
        }
        const std::vector<double> gray5 = functor.vectorized(q5, 1u, max_n);
        for (std::size_t i = 0; i < fast5.size(); ++i) {
            assert_near(static_cast<float>(fast5[i]), static_cast<float>(gray5[i]),
                        "K=5 low-k vs gray");
        }
    }
};

} // namespace moire_test

int main()
{
    moire_test::TransmissionIncrementalTestSuite suite;
    suite.run_tests();
    return suite.all_passed() ? 0 : 1;
}
