#pragma once
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"
#include "../VarRed/ControlVariate.h"
#include "MCPricer.h"

class VarRedMCPricer : public MCPricer
{
private:
    ControlVariate* CV;

public:
    VarRedMCPricer(RandomProcess* _process, PayOff* _payoff, ControlVariate* _cv, double _rate, double _maturity, size_t _nbSim, size_t _nbSteps);
    PriceResult Price() override;
};

