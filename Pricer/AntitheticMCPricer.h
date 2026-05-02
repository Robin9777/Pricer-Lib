#pragma once
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"
#include "../RandomGenerator/AntitheticNormal.h"
#include "MCPricer.h"

class AntitheticMCPricer : public MCPricer
{
    AntitheticNormal& antiNorm;

public:
    AntitheticMCPricer(RandomProcess* _process, PayOff* _payoff,
                       AntitheticNormal& _antiNorm,
                       double _rate, double _maturity,
                       size_t _nbSim, size_t _nbSteps);

    PriceResult Price() override;
};
