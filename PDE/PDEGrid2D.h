#pragma once
#include "R2R1Function.h"
#include "R1R1Function.h"
#include <vector>

class PDEGrid2D
{

private:
    void FillTopAndBottomBoundary();
    void FillRightBoundary();

protected:
    double T;
    double MinX;
    double MaxX;
    double h0;
    double h1;
    R2R1Function* a;
    R2R1Function* b;
    R2R1Function* r;
    R2R1Function* f;
    R1R1Function* TopBoundaryFunction;
    R1R1Function* BottomBoundaryFunction;
    R1R1Function* RightBoundaryFunction;
    std::vector<std::vector<double>> Nodes;
    size_t NodesHeight;
    size_t NodesWidth;

public:
    PDEGrid2D(
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

    virtual void FillNodes();
    double GetValue(double time, double spot);

};

