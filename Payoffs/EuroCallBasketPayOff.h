#pragma once
#include "PayOff.h"
class EuroCallBasketPayOff :
    public PayOff
{

private:
    double strike;
	std::vector<double> weights;

public:
    EuroCallBasketPayOff(double _strike, const std::vector<double>& _weights);
	double operator()(const std::vector<SinglePath*>& Paths) const override;

};

