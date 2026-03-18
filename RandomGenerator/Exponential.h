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
	UniformGenerator& Ugen;
	double GenerateInverseTransform() const;
	double GenerateRejectionSampling() const;


public:

	Exponential(double _lambda, UniformGenerator& _Ugen);
	virtual double Generate() override;
	double Generate(ExpoAlgo algo);

};

