#include "pch.h"
#include "PDEGridExplicit.h"


PDEGridExplicit::PDEGridExplicit(
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
    : PDEGrid2D(_T, _MinX, _MaxX, _NodesHeight, _NodesWidth,
        _a, _b, _r, _f,
        _TopBoundaryFunction, _BottomBoundaryFunction, _RightBoundaryFunction)
{
}

void PDEGridExplicit::FillNodes()
{

    for (size_t j = NodesWidth - 2; j != SIZE_MAX; j--) {

        for (size_t i = 1; i < NodesHeight - 1; i++) {

            double x = MinX + i * h0;
            double t = (j + 1) * h1;

            double A = (*a)(x, t);
            double B = (*b)(x, t);
            double R = (*r)(x, t);
            double F = (*f)(x, t);

            double Vup = Nodes[i + 1][j + 1];
            double Vmid = Nodes[i][j + 1];
            double Vdown = Nodes[i - 1][j + 1];

            double d2V = (Vup - 2.0 * Vmid + Vdown) / (h0 * h0);
            double dV = (Vup - Vdown) / (2.0 * h0);

            Nodes[i][j] = Vmid + h1 * (A * d2V + B * dV - R * Vmid - F);

        }

    }

}
