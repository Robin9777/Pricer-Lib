#pragma once
#include "BlackScholesND.h"

class BSMilsteinND : public BlackScholesND
{
public:
    BSMilsteinND(RandomGenerator* Gen,
                 const std::vector<double>& _spots,
                 double _rate,
                 const std::vector<double>& _vols,
                 const std::vector<std::vector<double>>& _corrMatrix);

    void Simulate(double startTime, double endTime, size_t nbSteps);
};
