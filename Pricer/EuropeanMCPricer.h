#pragma once
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"
#include "MCPricer.h"

class EuropeanMCPricer :
	public MCPricer
{
public:
	EuropeanMCPricer(
		RandomProcess* _process, 
		PayOff* _payoff,
		double _rate,
		double _maturity,
		size_t _nbSim,
		size_t _nbSteps
	);

	double Price() override;
};
