#pragma once
#include "PseudoGenerator.h"


class LinearCongruential :
    public PseudoGenerator
{

private:
    size_t multiplier;
    size_t increment;
    size_t modulus;
    // size_t seed  


public: 
    LinearCongruential() = default;
    LinearCongruential(size_t _seed, size_t _multiplier, size_t _increment, size_t _modulus);
    
    virtual double Generate();

    size_t GetModulus() const;
};

