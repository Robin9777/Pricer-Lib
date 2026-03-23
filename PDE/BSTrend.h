#pragma once
#include "R2R1Function.h"
class BSTrend :
    public R2R1Function
{

private:
    double rate;

public:
    BSTrend(double _rate);
    double operator()(double x, double t);
};

