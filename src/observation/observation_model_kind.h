#pragma once

#include <stdexcept>
#include <string>

enum class ObservationModelKind { Binary, Counts };

inline ObservationModelKind parse_observation_model_kind(const std::string &name)
{
    if (name == "binary")
    {
        return ObservationModelKind::Binary;
    }
    if (name == "counts")
    {
        return ObservationModelKind::Counts;
    }
    throw std::invalid_argument(
        "observation_model must be \"binary\" or \"counts\", not \"" + name + "\".");
}
