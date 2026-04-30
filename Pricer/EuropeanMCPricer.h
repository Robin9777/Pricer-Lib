#pragma once
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"

class EuropeanMCPricer
{

private:
	RandomProcess* Process; // BSMilstein1D for Vanilla option
	PayOff* Payoff;
	double maturity;
	double rate;
	size_t nbSim; // sim number
	size_t nbSteps; // steps number
	

public:
	EuropeanMCPricer(
		RandomProcess* process, 
		PayOff* payoff,
		double _maturity,
		double _rate,
		size_t _nbSim,
		size_t _nbSteps
	);

	double Price();

};

