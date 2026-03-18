#pragma once
#include "DiscreteGenerator.h"
#include <vector>
#include "UniformGenerator.h"


class FiniteSet :
    public DiscreteGenerator
{

private:

	std::vector<double> values;
	UniformGenerator& Ugen;

public:

	FiniteSet(std::vector<double> _values, UniformGenerator& _Ugen);

	virtual double Generate();

};

