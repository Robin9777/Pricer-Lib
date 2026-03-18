#pragma once
#define _USE_MATH_DEFINES
#include "pch.h"
#include "Normal.h"
#include <cmath>
#include "UniformGenerator.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Normal::Normal(double _mu, double _sigma, UniformGenerator& _Ugen) :
	mu(_mu),
	sigma(_sigma),
	Ugen(_Ugen)
{
}

double Normal::Generate() {
	return Generate(NormalAlgo::BoxMuller);
}

double Normal::Generate(NormalAlgo algo) {
	switch (algo) {
	case NormalAlgo::BoxMuller:
		return BoxMullerAlgorithm();
	case NormalAlgo::CentralLimitTheorem:
		return CentralLimitTheoremAlgorithm();
	case NormalAlgo::RejectionSampling:
		return RejectionSamplingAlgorithm();
	}
	return BoxMullerAlgorithm();
}

double Normal::BoxMullerAlgorithm() const {
	double U1 = Ugen.Generate();
	double U2 = Ugen.Generate();
	double Z0 = sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2);
	return this->mu + this->sigma * Z0;
}

double Normal::CentralLimitTheoremAlgorithm() const {
	double sum = 0.0;
	for (int i = 0; i < 12; ++i) {
		sum += Ugen.Generate();
	}
	return this->mu + this->sigma * (sum - 6.0);
}

double Normal::RejectionSamplingAlgorithm() const {
	double x, y;
	do {
		x = Ugen.Generate() * 10.0 - 5.0; // Generate x in range [-5, 5]
		y = Ugen.Generate() * (1.0 / (sigma * sqrt(2.0 * M_PI))); // Generate y in range [0, max density]
	} while (y > (1.0 / (sigma * sqrt(2.0 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2)));
	return x;
}