#pragma once
#include "DiscreteGenerator.h"
#include "UniformGenerator.h"


class Bernoulli :
    public DiscreteGenerator
{
private:
    UniformGenerator& Ugen;
    double p;

public:

    Bernoulli(const double& _p, UniformGenerator& _Ugen);

    // Methods
    virtual double Generate();


};

