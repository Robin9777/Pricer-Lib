#pragma once
#include "R1R1Function.h"
class CallTopBoundary :
    public R1R1Function
{

private:
    double Smax;
    double Strike;

public:

    CallTopBoundary(double _Smax, double _Strike);
    double operator()(double t) override;

};

