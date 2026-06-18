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
};
