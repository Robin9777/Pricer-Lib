#include "pch.h"
#include "EuropeanMCPricer.h"
#include <cmath>

EuropeanMCPricer::EuropeanMCPricer(RandomProcess* _process, PayOff* _payoff, double _rate, double _maturity, size_t _nbSim, size_t _nbSteps) :
	MCPricer(_process, _payoff, _rate, _maturity, _nbSim, _nbSteps)
{
}

double EuropeanMCPricer::Price()
{

	double totalPayoff = 0.0;
	double squaredPayoff = 0.0;


	for (size_t sim = 0; sim < nbSim; ++sim)
	{

		Process->Simulate(0.0, maturity, nbSteps);
		std::vector<SinglePath*> paths;

		for (int d = 0; d < Process->GetDimension(); ++d)
			paths.push_back(Process->GetPath(d));

		double discountedPayoff = std::exp(-rate * maturity) * (*Payoff)(paths);

		totalPayoff += discountedPayoff;
		squaredPayoff += discountedPayoff * discountedPayoff;
	}

	double averagePrice = totalPayoff / static_cast<double>(nbSim);
	double variance = (squaredPayoff / static_cast<double>(nbSim)) - (averagePrice * averagePrice);
	return averagePrice;
}
