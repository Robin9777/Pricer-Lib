#pragma once
#include "ContinuousGenerator.h"



enum class NormalAlgo {
    BoxMuller,
    CentralLimitTheorem,
    RejectionSampling
};


class Normal :
    public ContinuousGenerator
{

protected:
    // variables
    double mu;
	double sigma;

private:

	double BoxMullerAlgorythm() const;
	double CentralLimitTheoremAlgorythm() const;
	double RejectionSamplingAlgorythm() const;

public:

    Normal() = default;
    Normal(double _mu, double _sigma);

	virtual double Generate(NormalAlgo algo);



};

