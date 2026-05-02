#pragma once
#include "UniformGenerator.h"
#include <vector>

// Low-discrepancy generator based on the Halton / Van der Corput sequences.
//
// CORRECT USAGE FOR PATH SIMULATION:
//   For a T-step simulation using Box-Muller (2 uniforms per Normal draw):
//   - Provide T * normals_per_step * 2 prime bases (one per uniform needed)
//   - Use firstNPrimes(T * normals_per_step * 2) to generate them automatically
//
//   Example — 1D European, 252 steps, 1 Normal/step, Box-Muller:
//     HaltonGenerator h(HaltonGenerator::firstNPrimes(252 * 1 * 2));  // 504 primes
//
//   Example — 2D basket, 252 steps, 2 Normals/step:
//     HaltonGenerator h(HaltonGenerator::firstNPrimes(252 * 2 * 2));  // 1008 primes
//
// HOW IT WORKS:
//   The bases cycle: dim 0, 1, 2, ..., N-1, 0, 1, ...
//   Each full cycle of N dimensions = ONE simulation.
//   Simulation i uses Halton row i for ALL dimensions.
//   Different dimensions use different prime bases → approximately i.i.d. increments
//   within a path. Across simulations, each dimension's values follow the
//   low-discrepancy Van der Corput sequence → reduced variance.
//
// WHY NOT {2, 3} FOR MULTI-STEP?
//   With only 2 bases and T steps, the T Halton ROWS used within one simulation are
//   consecutive points of the same 2D sequence. These are anti-correlated by design
//   (low discrepancy in 2D), so sum(Z_j) ≈ 0 for every path → Var(W_T) ≈ 0.
//   The price is consistently wrong (biased). Using T*2 bases fixes this.
class HaltonGenerator : public UniformGenerator
{
    std::vector<size_t> bases;
    std::vector<size_t> counters;
    size_t currentDim;

    double radicalInverse(size_t n, size_t base) const;

public:
    explicit HaltonGenerator(size_t base = 2);
    explicit HaltonGenerator(const std::vector<size_t>& bases);

    double Generate() override;
    void Reset();

    // Returns the first n prime numbers, useful to build the bases vector.
    static std::vector<size_t> firstNPrimes(size_t n);
};
