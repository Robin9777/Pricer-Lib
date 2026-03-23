#include "pch.h"
#include "PDEGrid2D.h"


PDEGrid2D::PDEGrid2D(
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
)
    : T(_T),
    MinX(_MinX),
    MaxX(_MaxX),
    NodesHeight(_NodesHeight),
    NodesWidth(_NodesWidth),
    a(_a),
    b(_b),
    r(_r),
    f(_f),
    TopBoundaryFunction(_TopBoundaryFunction),
    BottomBoundaryFunction(_BottomBoundaryFunction),
    RightBoundaryFunction(_RightBoundaryFunction)
{

    h0 = (MaxX - MinX) / (NodesHeight - 1);
    h1 = T / (NodesWidth - 1);

    Nodes.assign(NodesHeight, std::vector<double>(NodesWidth, 0.0));
}


double PDEGrid2D::GetValue(double time, double spot)
{
    // To be implemented
    return 0.0;
}


void PDEGrid2D::FillTopAndBottomBoundary()
{
    for (size_t i = 0; i < NodesWidth; i++) {

        Nodes[0][i] = (*BottomBoundaryFunction)(i * h1);
        Nodes[NodesHeight - 1][i] = (*TopBoundaryFunction)(i * h1);

    }
}


void PDEGrid2D::FillRightBoundary()
{
    for (size_t i = 0; i < NodesHeight; i++) {
        
        double spot = MinX + i * h0;
        Nodes[i][NodesWidth - 1] = (*RightBoundaryFunction)(spot);

    }
} 