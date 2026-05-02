#pragma once
#include <cstddef>

class RandomGenerator
{

public:

	virtual double Generate() = 0;

	double Mean(size_t nbsim);
	double Variance(size_t nbsim);
};

