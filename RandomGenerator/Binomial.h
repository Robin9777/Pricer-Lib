#pragma once
#include "DiscreteGenerator.h"
#include "Bernoulli.h"
#include "UniformGenerator.h"

class Binomial :
    public DiscreteGenerator
{

private:
    Bernoulli bern;
    int n;


public:

    Binomial(int& _n, Bernoulli _bern);

    // Methods
    virtual double Generate();



};

