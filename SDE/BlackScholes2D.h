#pragma once
#include "RandomProcess.h"
#include "../RandomGenerator/RandomGenerator.h"
#include <vector>

class BlackScholes2D :
    public RandomProcess
{

protected:
    double Spot1;
    double Spot2;
    double Rate;
    double Vol1;
    double Vol2;
    std::vector<std::vector<double>> CorrelationMatrix;

public:
    BlackScholes2D(RandomGenerator* Gen,
                   double _spot1, double _spot2,
                   double _rate,
                   double _vol1, double _vol2,
                   double _rho);
};


