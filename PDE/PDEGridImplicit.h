#pragma once
#include "PDEGrid2D.h"
class PDEGridImplicit :
    public PDEGrid2D
{

public:
    PDEGridImplicit(
        double _T,
        double _MinX,
        double _MaxX,
        size_t _NodesHeight,
        size_t _NodesWidth,
        R2R1Function* _a,
        R2R1Function* _b,
        R2R1Function* _r,
        R2R1Function* _f,
        R1R1Function* _TopBoundaryFunction,
        R1R1Function* _BottomBoundaryFunction,
        R1R1Function* _RightBoundaryFunction
    );

    void FillNodes() override;
};

