#pragma once
#include "ContinuousGenerator.h"
#include "UniformGenerator.h"

enum class NormalAlgo {
	BoxMuller,
	CentralLimitTheorem,
	RejectionSampling
};

class Normal :
	public ContinuousGenerator
{
public:
	Normal() = delete;
	Normal(double _mu, double _sigma, UniformGenerator& _Ugen);

	virtual double Generate() override;
	double Generate(NormalAlgo algo);

protected:
	double mu;
	double sigma;
	UniformGenerator& Ugen;

private:
	double BoxMullerAlgorithm() const;
	double CentralLimitTheoremAlgorithm() const;
	double RejectionSamplingAlgorithm() const;
};
