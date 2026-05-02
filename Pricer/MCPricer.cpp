#include "pch.h"
#include "MCPricer.h"

MCPricer::MCPricer(RandomProcess* _process, PayOff* _payoff, double _rate, double _maturity, size_t _nbSim, size_t _nbSteps) :
	Process(_process),
	Payoff(_payoff),
	rate(_rate),
	maturity(_maturity),
	nbSim(_nbSim),
	nbSteps(_nbSteps)
{
}
