#include "ecoi_marginal.h"
#include "test_framework.hpp"

#include <cmath>
#include <span>
#include <vector>

namespace moire_test {

// Reference values are emitted by inst/scripts/ecoi_marginal_reference.R from the
// validated R oracle (inst/scripts/prototype_ecoi_marginal.R) for this exact
// deterministic scenario.
class EcoiMarginalTestSuite : public TestSuite {
public:
    EcoiMarginalTestSuite() : TestSuite("EcoiMarginal") {
        add_test("ecoi/r maps invert each other", [this]() { test_map_roundtrip(); });
        add_test("transmission term matches R oracle", [this]() { test_tx_matches_oracle(); });
        add_test("monoclonal atom is -inf for incompatible data",
                 [this]() { test_atom1(); });
        add_test("marginal L(e) matches R oracle", [this]() { test_marginal_matches_oracle(); });
        add_test("marginal LSE equals linear-space brute force",
                 [this]() { test_marginal_internal_consistency(); });
        add_test("conditional pi(m|e) matches oracle and normalizes",
                 [this]() { test_conditional_matches_oracle(); });
        add_test("population mixture matches oracle",
                 [this]() { test_population_mixture(); });
    }

private:
    static constexpr double kTol = 1e-3;

    static ecoi_marginal::Hyperparams hp()
    {
        return ecoi_marginal::Hyperparams{/*coi_p=*/0.5, /*coi_r=*/3.0,
                                          /*r_alpha=*/2.0, /*r_beta=*/3.0};
    }

    static constexpr int kMaxCoi = 8;

    static std::vector<ecoi_marginal::LocusSupport> scenario()
    {
        const std::vector<std::vector<float>> supports = {
            {0.2f, 0.3f},          // k = 2
            {0.1f, 0.15f, 0.25f},  // k = 3
            {0.4f},                // k = 1
        };
        std::vector<ecoi_marginal::LocusSupport> loci;
        for (const auto& s : supports) {
            loci.push_back(ecoi_marginal::precompute_locus(
                std::span<const float>(s), kMaxCoi));
        }
        return loci;
    }

    static void assert_near(double actual, double expected, const char* label)
    {
        if (std::fabs(actual - expected) > kTol) {
            throw std::runtime_error(
                std::string(label) + ": expected " + std::to_string(expected) +
                ", got " + std::to_string(actual));
        }
    }

    static void assert_neg_inf(double actual, const char* label)
    {
        if (!(std::isinf(actual) && actual < 0.0)) {
            throw std::runtime_error(std::string(label) + ": expected -inf, got " +
                                     std::to_string(actual));
        }
    }

    void test_map_roundtrip()
    {
        for (int m = 2; m <= 10; ++m) {
            for (double r : {0.05, 0.3, 0.5, 0.9}) {
                const double e = ecoi_marginal::ecoi_of(m, r);
                assert_near(ecoi_marginal::r_of(m, e), r, "r_of(ecoi_of)");
                // eCOI is bounded by COI and >= 1.
                ASSERT_GE(e, 1.0);
                ASSERT_LE(e, static_cast<double>(m));
            }
        }
    }

    void test_tx_matches_oracle()
    {
        const auto loci = scenario();
        assert_near(ecoi_marginal::tx_loglik_sample(loci, 3, 0.25),
                    -8.828889834896, "tx[m=3,r=0.25]");
        assert_near(ecoi_marginal::tx_loglik_sample(loci, 5, 0.60),
                    -8.975739553718, "tx[m=5,r=0.60]");
        // COI 2 cannot explain a locus with 3 distinct alleles.
        assert_neg_inf(ecoi_marginal::tx_loglik_sample(loci, 2, 0.10),
                       "tx[m=2,r=0.10]");
    }

    void test_atom1()
    {
        const auto loci = scenario();
        // Monoclonal is impossible here (locus 2 has 3 alleles).
        assert_neg_inf(ecoi_marginal::log_atom1(loci, hp()), "log_atom1");
    }

    void test_marginal_matches_oracle()
    {
        const auto loci = scenario();
        assert_near(ecoi_marginal::log_marginal_e(2.5, loci, hp(), kMaxCoi),
                    -10.238658017997, "log_marginal_e[2.5]");
        assert_near(ecoi_marginal::log_marginal_e(4.0, loci, hp(), kMaxCoi),
                    -12.662608448924, "log_marginal_e[4.0]");
        assert_near(ecoi_marginal::log_marginal_e(6.25, loci, hp(), kMaxCoi),
                    -18.372511602904, "log_marginal_e[6.25]");
        // e <= 1 has no continuous mass.
        assert_neg_inf(ecoi_marginal::log_marginal_e(1.0, loci, hp(), kMaxCoi),
                       "log_marginal_e[1.0]");
    }

    void test_marginal_internal_consistency()
    {
        const auto loci = scenario();
        const auto h = hp();
        for (double e : {2.5, 4.0, 6.25}) {
            const int m_lo = std::max(2, static_cast<int>(std::ceil(e)));
            long double linear = 0.0L;
            for (int m = m_lo; m <= kMaxCoi; ++m) {
                const double lw = ecoi_marginal::log_weight_m_given_e(m, e, loci, h);
                if (std::isfinite(lw)) {
                    linear += std::exp(static_cast<long double>(lw));
                }
            }
            const double brute = std::log(static_cast<double>(linear));
            assert_near(ecoi_marginal::log_marginal_e(e, loci, h, kMaxCoi), brute,
                        "marginal vs brute force");
        }
    }

    void test_conditional_matches_oracle()
    {
        const auto loci = scenario();
        const auto cond = ecoi_marginal::conditional_m(2.5, loci, hp(), kMaxCoi);

        const std::vector<int> exp_m = {3, 4, 5, 6, 7, 8};
        const std::vector<double> exp_r = {0.25, 0.5, 0.625, 0.7, 0.75, 0.785714285714};
        const std::vector<double> exp_w = {0.616993091885, 0.246686895379,
                                           0.086554297428, 0.032113132903,
                                           0.012549878077, 0.005102704328};
        ASSERT_EQ(exp_m.size(), cond.m.size());
        double wsum = 0.0;
        for (std::size_t i = 0; i < exp_m.size(); ++i) {
            ASSERT_EQ(exp_m[i], cond.m[i]);
            assert_near(cond.r[i], exp_r[i], "conditional r");
            assert_near(cond.w[i], exp_w[i], "conditional w");
            wsum += cond.w[i];
        }
        assert_near(wsum, 1.0, "conditional weights sum");

        // The conditional shares its terms with the marginal, so its
        // log-normalizer must equal log L(e).
        double approx_lse = 0.0;
        {
            const auto h = hp();
            std::vector<double> logw;
            for (int m : cond.m) {
                logw.push_back(ecoi_marginal::log_weight_m_given_e(m, 2.5, loci, h));
            }
            approx_lse =
                ecoi_marginal::log_sum_exp(std::span<const double>(logw));
        }
        assert_near(approx_lse,
                    ecoi_marginal::log_marginal_e(2.5, loci, hp(), kMaxCoi),
                    "conditional lse == marginal");
    }

    static std::vector<ecoi_marginal::LocusSupport> precompute_pop(
        const std::vector<std::vector<float>>& supports)
    {
        std::vector<ecoi_marginal::LocusSupport> loci;
        for (const auto& s : supports) {
            loci.push_back(ecoi_marginal::precompute_locus(
                std::span<const float>(s), kMaxCoi));
        }
        return loci;
    }

    void test_population_mixture()
    {
        const std::vector<double> log_pi = {std::log(0.6), std::log(0.4)};

        // Continuous part: same supports, different per-population frequencies.
        const std::vector<std::vector<ecoi_marginal::LocusSupport>> pops = {
            precompute_pop({{0.2f, 0.3f}, {0.1f, 0.15f, 0.25f}, {0.4f}}),
            precompute_pop({{0.3f, 0.3f}, {0.2f, 0.10f, 0.20f}, {0.5f}}),
        };
        assert_near(ecoi_marginal::sample_log_marginal_e(
                        std::span<const std::vector<ecoi_marginal::LocusSupport>>(pops),
                        std::span<const double>(log_pi), 2.5, hp(), kMaxCoi),
                    -9.695613900854, "sample_log_marginal_e[2.5]");
        assert_near(ecoi_marginal::sample_log_marginal_e(
                        std::span<const std::vector<ecoi_marginal::LocusSupport>>(pops),
                        std::span<const double>(log_pi), 4.0, hp(), kMaxCoi),
                    -11.825592789880, "sample_log_marginal_e[4.0]");

        // Monoclonal atom: all-k=1 loci so the atom is finite.
        const std::vector<std::vector<ecoi_marginal::LocusSupport>> pops_mono = {
            precompute_pop({{0.2f}, {0.3f}, {0.4f}}),
            precompute_pop({{0.3f}, {0.2f}, {0.5f}}),
        };
        assert_near(ecoi_marginal::sample_log_atom1(
                        std::span<const std::vector<ecoi_marginal::LocusSupport>>(pops_mono),
                        std::span<const double>(log_pi), hp()),
                    -5.174836309777, "sample_log_atom1_mono");
    }
};

}  // namespace moire_test

int main()
{
    moire_test::EcoiMarginalTestSuite suite;
    suite.run_tests();
    return suite.all_passed() ? 0 : 1;
}
