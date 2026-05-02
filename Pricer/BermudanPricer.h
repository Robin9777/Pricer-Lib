#pragma once
#include "MCPricer.h"
#include "../Payoffs/Payoff.h"
#include "../SDE/RandomProcess.h"
#include "Eigen/Dense"
#include <vector>

class BermudanPricer :
    public MCPricer
{

protected:
    // Berm Members
    std::vector<double> exerciseDates;
    const int PolynomialDegree = 3;

private:
    std::vector<std::vector<int>> generateBasisDegrees(int D, int p);
    Eigen::VectorXd regressionFit(const std::vector<std::vector<double>>& X, const std::vector<double>& Y, const std::vector<std::vector<int>>& basisDegrees);

public:
    BermudanPricer(RandomProcess* _process, PayOff* _payoff, double _rate, double _maturity, size_t _nbSim, size_t _nbSteps, const std::vector<double>& _exerciseDates);
    double Price() override;
    virtual ~BermudanPricer();

};
