#include "pch.h"
#include "BSMilsteinND.h"
#include "BrownianND.h"
#include "../SDE/SinglePath.h"
#include <cmath>

BSMilsteinND::BSMilsteinND(RandomGenerator* Gen,
    const std::vector<double>& _spots,
    double _rate,
    const std::vector<double>& _vols,
    const std::vector<std::vector<double>>& _corrMatrix)
    : BlackScholesND(Gen, _spots, _rate, _vols, _corrMatrix)
{
}

void BSMilsteinND::Simulate(double startTime, double endTime, size_t nbSteps)
{
    int N = this->Dimension;

    for (auto p : this->Paths) delete p;
    this->Paths.clear();
    for (int d = 0; d < N; d++)
        this->Paths.push_back(new SinglePath(startTime, endTime, nbSteps));

    BrownianND W(this->Generator, N, &this->CorrelationMatrix);
    W.Simulate(startTime, endTime, nbSteps);

    double dt = (endTime - startTime) / nbSteps;

    std::vector<double> S(this->Spots);
    for (int d = 0; d < N; d++)
        this->Paths[d]->InsertValue(S[d]);

    for (size_t i = 1; i <= nbSteps; i++) {
        double t_prev = startTime + (i - 1) * dt;
        double t_curr = startTime + i * dt;

        for (int d = 0; d < N; d++) {
            double dW = W.GetPath(d)->GetState(t_curr) - W.GetPath(d)->GetState(t_prev);
            S[d] += this->Rate * S[d] * dt
                  + this->Vols[d] * S[d] * dW
                  + 0.5 * std::pow(this->Vols[d], 2) * S[d] * (std::pow(dW, 2) - dt);
            this->Paths[d]->InsertValue(S[d]);
        }
    }
}
