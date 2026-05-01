#pragma once
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"

class EuropeanMCPricer
{

private:
	RandomProcess* Process; // BSMilstein1D for Vanilla option
	PayOffTemplate* POTemplate;
	PayOff* Payoff; // EuropeanCallPayoff for Vanilla option
	// double rate; already defined in  BSMilstein
	size_t nbSim; // sim number
	size_t nbSteps; // steps number
	

public:
	EuropeanMCPricer(
		RandomProcess* _process, 
		PayOffTemplate* _payofftemplate,
		PayOff* _payoff,
		size_t _nbSim,
		size_t _nbSteps
	);

	double Price();

};

