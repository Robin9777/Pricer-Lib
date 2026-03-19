#pragma once
#include "VanillaTerminalCondition.h"
class CallTerminalCondition :
    public VanillaTerminalCondition
{

public:
    CallTerminalCondition() = default;
    CallTerminalCondition(double _strike);

    double operator()(double x) override;
    
};

