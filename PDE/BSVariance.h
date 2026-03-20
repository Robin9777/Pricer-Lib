#pragma once
#include "R2R1Function.h"
class BSVariance :
    public R2R1Function
{
private:
    double Sigma;

public:
    BSVariance(double _sigma);
    double operator()(double x, double t) override;
};

