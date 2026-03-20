#pragma once
#include "R1R1Function.h"
class CallBottomBoundary :
    public R1R1Function
{

private:
    double Smin;
    double Strike;

public:
    CallBottomBoundary(double _Smin, double _Strike);
    double operator()(double t) override;
};

