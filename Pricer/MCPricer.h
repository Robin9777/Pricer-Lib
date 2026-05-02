#pragma once
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"

class MCPricer
{
protected:
    RandomProcess* Process;
    PayOff* Payoff;
    double         rate;
    double         maturity;
    size_t         nbSim;
    size_t         nbSteps;

public:
    MCPricer(RandomProcess* _process, PayOff* _payoff, double _rate, double _maturity, size_t _nbSim, size_t _nbSteps);
    virtual double Price() = 0;
    virtual ~MCPricer() = default;
};

