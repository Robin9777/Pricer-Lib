#pragma once
#include "VanillaTerminalCondition.h"
class PutTerminalCondition :
    public VanillaTerminalCondition
{
public:
    PutTerminalCondition(double _strike);
    double operator()(double x) override;
};

