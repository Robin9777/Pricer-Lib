#pragma once
#include "BlackScholes1D.h"
class BSEuler1D :
    public BlackScholes1D
{

public:
    BSEuler1D(RandomGenerator* Gen, double _spot, double _rate, double _vol);
    void Simulate(double startTime, double endTime, size_t nbSteps);
};

