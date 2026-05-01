#pragma once
#include "PayOff.h"
class EuropeanCallPayoff :
    public PayOff
{

private:
    double strike;

public:
    EuropeanCallPayoff(double _strike);
    double operator()(const std::vector<SinglePath*>& Paths) const override;


};

