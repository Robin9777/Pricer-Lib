#pragma once
#include "R1R1Function.h"
class VanillaTerminalCondition :
    public R1R1Function
{

protected:
    double Strike;

public:
    VanillaTerminalCondition() = default;
    VanillaTerminalCondition(double _strike);

};

