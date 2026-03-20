#pragma once
#include "R2R1Function.h"
class NullFunction :
    public R2R1Function
{

public:
    NullFunction();
    double operator()(double x, double t) override;
};