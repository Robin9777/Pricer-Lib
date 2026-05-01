#pragma once
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"

class EuropeanMCPricer
{

private:
	RandomProcess* Process; // BSMilstein1D for Vanilla option
	PayOff* Payoff; // EuropeanCallPayoff for Vanilla option
	double rate;
	double maturity;
	size_t nbSim; // sim number
	size_t nbSteps; // steps number
	

public:
	EuropeanMCPricer(
		RandomProcess* _process, 
		PayOff* _payoff,
		double _rate,
		double _maturity,
		size_t _nbSim,
		size_t _nbSteps
	);

	double Price();

};