#pragma once
#include "RandomProcess.h"
#include <vector>
class BrownianND :
    public RandomProcess
{

protected:
    std::vector<std::vector<double>>* CorrelationMatrix;
    std::vector<std::vector<double>> cholesky(const std::vector<std::vector<double>>& matrix);

public:
    BrownianND(RandomGenerator* Gen, int dim, std::vector<std::vector<double>>* Corr);
    void Simulate(double startTime, double endTime, size_t nbSteps);

};
