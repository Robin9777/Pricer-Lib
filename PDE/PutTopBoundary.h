#pragma once
#include "R1R1Function.h"
class PutTopBoundary :
    public R1R1Function
{

private:
    double Smax;
    double Strike;

public:

    PutTopBoundary(double _Smax, double _Strike);
    double operator()(double t) override;

};