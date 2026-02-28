#pragma once
#include "DiscreteGenerator.h"
#include <vector>
#include "Bernoulli.h"


class FiniteSet :
    public DiscreteGenerator
{

private:

	std::vector<double> values;
	Bernoulli Bern;

public:

	FiniteSet() = default;
	FiniteSet(std::vector<double> _values, Bernoulli _Bern);

	virtual double Generate();

};

