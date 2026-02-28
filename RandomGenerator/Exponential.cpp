#pragma once
#include "Exponential.h"
#include <cmath>
#include "UniformGenerator.h"


Exponential::Exponential(double _lambda, UniformGenerator& _Ugen) :
	lambda(_lambda)
{
}

Exponential::Generate(ExpoAlgo algo) {
	switch (algo) {
		case ExpoAlgo::InverseTransform:
			return GenerateInverseTransform();
		case ExpoAlgo::RejectionSampling:
			return GenerateRejectionSampling();
	}
}

Exponential::GenerateInverseTransform() const {
	double U = UniformGenerator::Generate();
	return -log(1.0 - U) / lambda;
}

Exponential::GenerateRejectionSampling() const {
	double x, y;
	do {
		x = UniformGenerator::Generate() * 10.0; // Generate x in range [0, 10]
		y = UniformGenerator::Generate() * lambda; // Generate y in range [0, lambda]
	} while (y > lambda * exp(-lambda * x));
	return x;
}