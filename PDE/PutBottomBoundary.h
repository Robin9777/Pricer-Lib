#pragma once
#include "R1R1Function.h"
class PutBottomBoundary :
    public R1R1Function
{
private:
    double Smin;
    double Strike;

public:
    PutBottomBoundary(double _Smin, double _Strike);
    double operator()(double t) override;
};

