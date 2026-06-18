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
