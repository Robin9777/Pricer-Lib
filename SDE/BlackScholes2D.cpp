#include "pch.h"
#include "BlackScholes2D.h"

BlackScholes2D::BlackScholes2D(RandomGenerator* Gen,
    double _spot1, double _spot2,
    double _rate,
    double _vol1, double _vol2,
    double _rho) :
    RandomProcess(Gen, 2),
    Spot1(_spot1),
    Spot2(_spot2),
    Rate(_rate),
    Vol1(_vol1),
    Vol2(_vol2),
    CorrelationMatrix({ {1.0, _rho}, {_rho, 1.0} })
{
}
