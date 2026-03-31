#pragma once
#include "RandomProcess.h"
#include "../RandomGenerator/RandomGenerator.h"
#include <vector>
class Heston :
    public RandomProcess
{

protected:
    double Spot;
    double InitVariance;
    double Mu;
    double Theta;
    double Kappa;
    double Sigma;
    double Rho;
	std::vector<std::vector<double>> CorrelationMatrix;

public:
    Heston(RandomGenerator* Gen,
           double _spot,
           double _initVariance,
           double _mu,
           double _theta,
           double _kappa,
           double _sigma,
		   double _rho
    );

	void Simulate(double startTime, double endTime, size_t nbSteps);
};

