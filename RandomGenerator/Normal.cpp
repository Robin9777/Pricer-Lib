#pragma once
#include "Normal.h"
#include <cmath>
#include "UniformGenerator.h"



Normal::Normal(double _mu, double _sigma) :
	mu(_mu),
	sigma(_sigma)
{

}

double Normal::Generate(NormalAlgo algo) {
	switch (algo) {
		case NormalAlgo::BoxMuller:
			return BoxMullerAlgorythm();
		case NormalAlgo::CentralLimitTheorem:
			return CentralLimitTheoremAlgorythm();
		case NormalAlgo::RejectionSampling:
			return RejectionSamplingAlgorythm();
	}
}


double Normal::BoxMullerAlgorythm() const {
	double U1 = UniformGenerator::Generate();
	double U2 = UniformGenerator::Generate();
	double Z0 = sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2);
	return this->mu + this->sigma * Z0;
}

double Normal::CentralLimitTheoremAlgorythm() const {

	// Why 12? Why not 10000?
	double sum = 0.0;
	for (int i = 0; i < 12; ++i) {
		sum += UniformGenerator::Generate();
	}
	return this->mu + this->sigma * (sum - 6.0);
}

double Normal::RejectionSamplingAlgorythm() const {
	double x, y;
	do {
		x = UniformGenerator::Generate() * 10.0 - 5.0; // Generate x in range [-5, 5]
		y = UniformGenerator::Generate() * (1.0 / (sigma * sqrt(2.0 * M_PI))); // Generate y in range [0, max density]
	} while (y > (1.0 / (sigma * sqrt(2.0 * M_PI))) * exp(-0.5 * pow((x - mu) / sigma, 2)));
	return x;
}