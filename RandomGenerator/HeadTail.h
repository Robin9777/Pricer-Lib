#pragma once
#include "DiscreteGenerator.h"
#include "UniformGenerator.h"

class HeadTail :
    public DiscreteGenerator
{

private:
    UniformGenerator& Ugen;

public:

    HeadTail() = default;
    HeadTail(UniformGenerator& _Ugen);

    // Methods
    virtual double Generate();

};

