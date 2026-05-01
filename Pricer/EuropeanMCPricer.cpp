#include "pch.h"
#include "EuropeanMCPricer.h"
#include <cmath>

EuropeanMCPricer::EuropeanMCPricer(RandomProcess* _process, PayOffTemplate* _payofftemplate, PayOff* _payoff, size_t _nbSim, size_t _nbSteps) :
	Process(_process),
	POTemplate(_payofftemplate),
	Payoff(_payoff),
	nbSim(_nbSim),
	nbSteps(_nbSteps)
{
}

double EuropeanMCPricer::Price()
{
	const double maturity = POTemplate->maturity;
	const double rate     = Process->GetRate();

	double totalPayoff = 0.0;

	for (size_t sim = 0; sim < nbSim; ++sim)
	{
		POTemplate->Paths.clear();

		Process->Simulate(0.0, maturity, nbSteps);

		for (int d = 0; d < Process->GetDimension(); ++d)
			POTemplate->Paths.push_back(Process->GetPath(d));


		totalPayoff += (*Payoff)(*POTemplate);
	}

	double avgPayoff = totalPayoff / static_cast<double>(nbSim);
	return std::exp(-rate * maturity) * avgPayoff;
}
