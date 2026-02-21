#pragma once
#include "UniformGenerator.h"

class PseudoGenerator : public UniformGenerator
{
    
private:
    size_t seed;

public:

    PseudoGenerator() = default;
    PseudoGenerator(size_t _seed);

    // getter
    size_t GetSeed() const {
        return seed;
    }

    // setter
    void SetSeed(const size_t& seed_) {
        this->seed = seed_;
    }
};

