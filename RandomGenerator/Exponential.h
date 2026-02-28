#pragma once
#include "ContinuousGenerator.h"
#include "UniformGenerator.h"


enum class ExpoAlgo {
	InverseTransform,
	RejectionSampling
};


class Exponential :
    public ContinuousGenerator
{

private:

	double lambda;
	double GenerateInverseTransform() const;
	double GenerateRejectionSampling() const;


public:

	Exponential() = default;
	Exponential(double _lambda, UniformGenerator& _Ugen);
	virtual double Generate(ExpoAlgo algo);

};

