#pragma once
#include "RandomProcess.h"
#include "../RandomGenerator/RandomGenerator.h"
class BlackScholes1D :
    public RandomProcess
{

protected:
    double Spot;
    double Rate;
    double Vol;

public:
    BlackScholes1D(RandomGenerator* Gen, double _spot, double _rate, double _vol);
};

