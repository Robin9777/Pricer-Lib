#pragma once
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"
#include "../VarRed/ControlVariate.h"

class VarRedMCPricer
{

private:
    RandomProcess* Process;
    PayOff* Payoff;
    ControlVariate* CV;
    double Rate;
    double Maturity;
    size_t NbSim;
    size_t NbSteps;

public:
    VarRedMCPricer(RandomProcess*, PayOff*, ControlVariate*,
        double rate, double maturity, size_t nbSim, size_t nbSteps);
    double Price();
};

