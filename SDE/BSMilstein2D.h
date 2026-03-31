#pragma once
#include "BlackScholes2D.h"

class BSMilstein2D :
    public BlackScholes2D
{

public:
    BSMilstein2D(RandomGenerator* Gen,
                 double _spot1, double _spot2,
                 double _rate,
                 double _vol1, double _vol2,
                 double _rho);
    void Simulate(double startTime, double endTime, size_t nbSteps);
};


