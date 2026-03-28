#pragma once
#include "RandomProcess.h"
#include "../RandomGenerator/RandomGenerator.h"
class Brownian1D :
    public RandomProcess
{

public:
    Brownian1D(RandomGenerator* Gen);
    void Simulate(double startTime, double endTime, size_t nbSteps);

};

