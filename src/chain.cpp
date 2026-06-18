#include "chain.h"

#include "prob_any_missing.h"
#include "multivector_fused.h"
#include "mcmc_utils.h"
#include "sampler.h"
#include "profiler.h"

#include <cmath>
#include <algorithm>

#include "env_defs.h"

#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <span>

constexpr float min_sampled = std::numeric_limits<float>::min();

void Chain::set_llik(float llik) { this->llik = llik; }

void Chain::set_temp(float temp) { this->temp = temp; }

float Chain::get_temp() { return this->temp; }

Chain::Chain(GenotypingData genotyping_data, Parameters params, float temp)
    : genotyping_data(genotyping_data),
      params(params),
      observation_model_(make_observation_model(params.observation_model_kind,
                                                params.max_eps_neg,
                                                params.max_eps_pos)),
      sampler(),
      temp(temp),
      llik(std::numeric_limits<float>::lowest())

{
    genotyping_data.validate_for_observation_model(params.observation_model_kind);
    initialize_parameters();
}
