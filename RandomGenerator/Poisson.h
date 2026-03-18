#pragma once
#include "DiscreteGenerator.h"
#include "UniformGenerator.h"


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
	UniformGenerator& Ugen;

    double GenerateAlgo1();
    double GenerateAlgo2();

public:

    Poisson(double _lambda, UniformGenerator& _Ugen);

    virtual double Generate() override;
    double Generate(PoissonAlgo algo);

};

