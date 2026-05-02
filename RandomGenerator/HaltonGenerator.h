#pragma once
#include "UniformGenerator.h"
#include <vector>

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

    static std::vector<size_t> firstNPrimes(size_t n);
};
