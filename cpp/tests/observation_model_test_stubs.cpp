#include "binary_observation_model.h"
#include "count_poisson_observation_model.h"
#include "sampler.h"

// Test build links log-likelihood only; sampling is exercised via R integration tests.
LatentGenotype BinaryObservationModel::sample_latent_genotype(
    Sampler&, std::span<int const>, int, float, float) const
{
    return LatentGenotype{{0}, 0.0f};
}

LatentGenotype CountPoissonObservationModel::sample_latent_genotype(
    Sampler&, std::span<int const>, int, float, float) const
{
    return LatentGenotype{{0}, 0.0f};
}

// Marginal-proposal methods are exercised at the observation_model_math layer in
// the unit tests; these stubs keep the test link self-contained (no boost).
void BinaryObservationModel::marginal_presence_probs(std::span<int const>, float, float,
                                                     std::vector<float>&) const
{
}

LatentGenotype BinaryObservationModel::propose_latent_genotype_marginal(
    Sampler&, std::span<int const>, float, float) const
{
    return LatentGenotype{{0}, 0.0f};
}

float BinaryObservationModel::latent_genotype_log_prob_marginal(
    std::span<int const>, std::span<int const>, float, float) const
{
    return 0.0f;
}

void CountPoissonObservationModel::marginal_presence_probs(std::span<int const>, float, float,
                                                           std::vector<float>&) const
{
}

LatentGenotype CountPoissonObservationModel::propose_latent_genotype_marginal(
    Sampler&, std::span<int const>, float, float) const
{
    return LatentGenotype{{0}, 0.0f};
}

float CountPoissonObservationModel::latent_genotype_log_prob_marginal(
    std::span<int const>, std::span<int const>, float, float) const
{
    return 0.0f;
}
