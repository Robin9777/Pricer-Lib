#include "pch.h"
#include "BlackScholesND.h"

BlackScholesND::BlackScholesND(RandomGenerator* Gen,
    const std::vector<double>& _spots,
    double _rate,
    const std::vector<double>& _vols,
    const std::vector<std::vector<double>>& _corrMatrix)
    : RandomProcess(Gen, (int)_spots.size()),
      Spots(_spots),
      Rate(_rate),
      Vols(_vols),
      CorrelationMatrix(_corrMatrix)
{
}
