#pragma once

#include "latent_genotype.h"

#include <span>

class Sampler;

class ObservationModel {
   public:
    virtual ~ObservationModel() = default;

    virtual float log_likelihood(std::span<int const> latent_allele_indices,
                                 std::span<int const> observed_barcode,
                                 float epsilon_neg,
                                 float epsilon_pos) const = 0;

    virtual LatentGenotype sample_latent_genotype(Sampler &sampler,
                                                  std::span<int const> observed_barcode,
                                                  int coi,
                                                  float epsilon_pos,
                                                  float epsilon_neg) const = 0;

    // COI-independent latent-genotype proposal for the fully-marginalized (eCOI)
    // sampler. Returns a non-empty genotype and the log proposal density
    // log q(G | observed, eps); the density is exactly reproducible by
    // latent_genotype_log_prob_marginal so the reverse move can be scored.
    virtual LatentGenotype propose_latent_genotype_marginal(
        Sampler &sampler, std::span<int const> observed_barcode, float epsilon_pos,
        float epsilon_neg) const = 0;

    // Evaluate log q(G | observed, eps) of the marginal proposal for an arbitrary
    // (possibly -1 padded) latent genotype. Must match the density returned by
    // propose_latent_genotype_marginal.
    virtual float latent_genotype_log_prob_marginal(
        std::span<int const> latent_allele_indices, std::span<int const> observed_barcode,
        float epsilon_pos, float epsilon_neg) const = 0;
};
