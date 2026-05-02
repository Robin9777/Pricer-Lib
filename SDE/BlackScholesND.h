#pragma once
#include "RandomProcess.h"
#include "../RandomGenerator/RandomGenerator.h"
#include <vector>

class BlackScholesND : public RandomProcess
{
protected:
    std::vector<double> Spots;
    double Rate;
    std::vector<double> Vols;
    std::vector<std::vector<double>> CorrelationMatrix;

public:
    BlackScholesND(RandomGenerator* Gen,
                   const std::vector<double>& _spots,
                   double _rate,
                   const std::vector<double>& _vols,
                   const std::vector<std::vector<double>>& _corrMatrix);

    double GetRate() const override { return Rate; }
};
