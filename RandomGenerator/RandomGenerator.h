#pragma once

class RandomGenerator
{

public:

	virtual double Generate() = 0;

	double Mean(size_t nbsim);
	double Variance(size_t nbsim);
};

