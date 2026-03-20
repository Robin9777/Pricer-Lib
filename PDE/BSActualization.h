#pragma once
#include "R2R1Function.h"
class BSActualization :
    public R2R1Function
{
private:
    double Rate;

public:
    BSActualization(double _Rate);
    double operator()(double x, double t) override;
};

