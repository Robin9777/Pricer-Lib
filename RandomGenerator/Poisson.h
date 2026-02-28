#pragma once
#include "DiscreteGenerator.h"
#include "Bernoulli.h"


enum class PoissonAlgo {
    PoissonAlgo1,
    PoissonAlgo2
};

class Poisson :
    public DiscreteGenerator
{

private:

    // variables
    double lambda;
	Bernoulli Bern;

    double GenerateAlgo1();
    double GenerateAlgo2();

public:

    Poisson() = default;
    Poisson(double _lambda, Bernoulli _Bern);

    virtual double Generate(PoissonAlgo algo);

};

