#include "pch.h"
#include "EuropeanMCPricer.h"

EuropeanMCPricer::EuropeanMCPricer(RandomProcess* _process, PayOff* _payoff, double _maturity, double _rate, size_t _nbSim, size_t _nbSteps)
	: Process(_process),
	Payoff(_payoff),
	maturity(_maturity),
	rate(_rate),
	nbSim(_nbSim),
	nbSteps(_nbSteps)
{
}

double EuropeanMCPricer::Price()
{


	return 0.0;
}
