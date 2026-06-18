#pragma once

#include "binary_observation_model.h"
#include "count_poisson_observation_model.h"
#include "observation_model.h"
#include "observation_model_kind.h"

#include <memory>
#include <stdexcept>

inline std::unique_ptr<ObservationModel> make_observation_model(ObservationModelKind kind,
                                                                float max_eps_neg,
                                                                float max_eps_pos)
{
    switch (kind)
    {
        case ObservationModelKind::Binary:
            return std::make_unique<BinaryObservationModel>(max_eps_neg, max_eps_pos);
        case ObservationModelKind::Counts:
            return std::make_unique<CountPoissonObservationModel>(max_eps_neg, max_eps_pos);
    }
    throw std::logic_error("Unhandled observation model kind.");
}
