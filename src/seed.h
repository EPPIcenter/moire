#pragma once

#ifndef SEED_H_
#define SEED_H_

#include <cstddef>
#include <cstdint>

// Child RNG derivation for run_mcmc(). Independent MCMC chains (R num_chains)
// are spaced independent_chain_stride apart so PT replica offsets (parent + 0,
// 1, ...) cannot collide. Keep independent_chain_stride in sync with
// resolve_mcmc_seeds() in R/mcmc.R.
namespace Seed
{
inline constexpr std::uint32_t independent_chain_stride = 1000;
inline constexpr std::uint32_t swap_rng_salt = 0x9E3779B9u;

inline std::uint32_t pt_replica(std::uint32_t parent, std::size_t replica)
{
    return parent + static_cast<std::uint32_t>(replica);
}

inline std::uint32_t swap_sampler(std::uint32_t parent)
{
    return parent ^ swap_rng_salt;
}
}  // namespace Seed

#endif  // SEED_H_
