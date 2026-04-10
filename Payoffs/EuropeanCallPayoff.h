#pragma once
#include "PayOff.h"
class EuropeanCallPayoff :
    public PayOff
{
public:
    double operator()(PayOffTemplate PO_template) override;
};

