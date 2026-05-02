#include "pch.h"
#include "AntitheticMCPricer.h"
#include <cmath>

AntitheticMCPricer::AntitheticMCPricer(RandomProcess* _process, PayOff* _payoff,
                                        AntitheticNormal& _antiNorm,
                                        double _rate, double _maturity,
                                        size_t _nbSim, size_t _nbSteps)
    : MCPricer(_process, _payoff, _rate, _maturity, _nbSim, _nbSteps),
      antiNorm(_antiNorm)
{
}

PriceResult AntitheticMCPricer::Price()
{
    double totalPayoff   = 0.0;
    double squaredPayoff = 0.0;
    double discount      = std::exp(-rate * maturity);

    for (size_t sim = 0; sim < nbSim; ++sim) {

        // --- Normal path ---
        antiNorm.ResetBuffer();
        antiNorm.SetAntithetic(false);
        Process->Simulate(0.0, maturity, nbSteps);
        std::vector<SinglePath*> paths;
        for (int d = 0; d < Process->GetDimension(); ++d)
            paths.push_back(Process->GetPath(d));
        double y = discount * (*Payoff)(paths);

        // --- Antithetic path (same Z values, negated) ---
        antiNorm.ResetIndex();
        antiNorm.SetAntithetic(true);
        Process->Simulate(0.0, maturity, nbSteps);
        std::vector<SinglePath*> paths_anti;
        for (int d = 0; d < Process->GetDimension(); ++d)
            paths_anti.push_back(Process->GetPath(d));
        double y_anti = discount * (*Payoff)(paths_anti);

        double avg = 0.5 * (y + y_anti);
        totalPayoff   += avg;
        squaredPayoff += avg * avg;
    }

    double mean     = totalPayoff / static_cast<double>(nbSim);
    double variance = (squaredPayoff / static_cast<double>(nbSim)) - mean * mean;
    return { mean, 1.96 * std::sqrt(variance / static_cast<double>(nbSim)) };
}
