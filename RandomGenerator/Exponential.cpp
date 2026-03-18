#pragma once
#include "pch.h"
#include "Exponential.h"
#include <cmath>
#include "UniformGenerator.h"


Exponential::Exponential(double _lambda, UniformGenerator& _Ugen) :
	lambda(_lambda),
	Ugen(_Ugen)
{
}

double Exponential::Generate() {
	return Generate(ExpoAlgo::InverseTransform);
}

double Exponential::Generate(ExpoAlgo algo) {
	switch (algo) {
	case ExpoAlgo::InverseTransform:
		return GenerateInverseTransform();
	case ExpoAlgo::RejectionSampling:
		return GenerateRejectionSampling();
	}
}

double Exponential::GenerateInverseTransform() const {
	double U = Ugen.Generate();
	return -log(1.0 - U) / lambda;
}

double Exponential::GenerateRejectionSampling() const {
	double x, y;
	do {
		x = Ugen.Generate() * 10.0; // Generate x in range [0, 10]
		y = Ugen.Generate() * lambda; // Generate y in range [0, lambda]
	} while (y > lambda * exp(-lambda * x));
	return x;
}