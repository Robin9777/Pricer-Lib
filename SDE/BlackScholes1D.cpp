#include "pch.h"
#include "BlackScholes1D.h"

BlackScholes1D::BlackScholes1D(RandomGenerator* Gen, double _spot, double _rate, double _vol) : 
	RandomProcess(Gen, 1), 
	Spot(_spot),
	Rate(_rate),
	Vol(_vol)
{
}
