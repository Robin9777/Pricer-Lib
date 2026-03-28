#include "pch.h"
#include "PDEGridImplicit.h"


PDEGridImplicit::PDEGridImplicit(
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

void PDEGridImplicit::FillNodes()
{
    size_t M = NodesHeight;  // space
    size_t N = NodesWidth;   // time

    double dx = h0;
    double dt = h1;

    // Vectors for Thomas algorithm
    std::vector<double> lower(M);
    std::vector<double> diag(M);
    std::vector<double> upper(M);
    std::vector<double> rhs(M);

    // Backward in time
    for (size_t j = N - 2; j != SIZE_MAX; j--)
    {
        double t = (j + 1) * dt;

        // Build system
        for (size_t i = 1; i < M - 1; i++)
        {
            double x = MinX + i * dx;

            double A = (*a)(x, t);
            double B = (*b)(x, t);
            double R = (*r)(x, t);
            double F = (*f)(x, t);

            double alpha =
                -dt * (A / (dx * dx)
                    - B / (2.0 * dx));

            double beta =
                1.0
                + dt * (2.0 * A / (dx * dx)
                    + R);

            double gamma =
                -dt * (A / (dx * dx)
                    + B / (2.0 * dx));

            lower[i] = alpha;
            diag[i] = beta;
            upper[i] = gamma;

            rhs[i] =
                Nodes[i][j + 1]
                + dt * F;
        }

        // Boundary conditions
        rhs[1] -= lower[1] * Nodes[0][j];
        rhs[M - 2] -= upper[M - 2] * Nodes[M - 1][j];

        // Thomas algorithm

        // Forward elimination
        for (size_t i = 2; i < M - 1; i++)
        {
            double m =
                lower[i] / diag[i - 1];

            diag[i] -= m * upper[i - 1];

            rhs[i] -= m * rhs[i - 1];
        }

        // Back substitution
        Nodes[M - 2][j] =
            rhs[M - 2] / diag[M - 2];

        for (size_t i = M - 3;
            i != SIZE_MAX;
            i--)
        {
            if (i == 0)
                break;

            Nodes[i][j] =
                (rhs[i]
                    - upper[i]
                    * Nodes[i + 1][j])
                / diag[i];
        }
    }
}