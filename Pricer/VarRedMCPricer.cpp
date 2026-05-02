#include "pch.h"
#include "VarRedMCPricer.h"
#include <cmath>


VarRedMCPricer::VarRedMCPricer(RandomProcess* _process, PayOff* _payoff, ControlVariate* _cv, double _rate, double _maturity, size_t _nbSim, size_t _nbSteps) :
	MCPricer(_process, _payoff, _rate, _maturity, _nbSim, _nbSteps), CV(_cv)
{
}

PriceResult VarRedMCPricer::Price()
{
	double totalPayoff = 0.0;
	double squaredPayoff = 0.0;

	double cvExpectation = CV->AnalyticalExpectation();

	for (size_t sim = 0; sim < nbSim; ++sim) {
		
		//classic approach
		Process->Simulate(0.0, maturity, nbSteps);
		std::vector<SinglePath*> paths;
		for (int d = 0; d < Process->GetDimension(); ++d)
			paths.push_back(Process->GetPath(d));

		double discountedPayoff = std::exp(-rate * maturity) * (*Payoff)(paths);
		double cvValue = std::exp(-rate * maturity) * CV->SimulatedValue(paths);
		double adjustedPayoff = discountedPayoff - cvValue + cvExpectation;

		totalPayoff += adjustedPayoff;
		squaredPayoff += adjustedPayoff * adjustedPayoff;
	}

	double averagePrice = totalPayoff / static_cast<double>(nbSim);
	double variance = (squaredPayoff / static_cast<double>(nbSim)) - (averagePrice * averagePrice);

	return {averagePrice, 1.96 * std::sqrt(variance / static_cast<double>(nbSim))};
}
