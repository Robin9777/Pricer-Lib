#include "pch.h"
#include "BSMilstein2D.h"
#include "BrownianND.h"
#include "../SDE/SinglePath.h"
#include <cmath>

BSMilstein2D::BSMilstein2D(RandomGenerator* Gen,
    double _spot1, double _spot2,
    double _rate,
    double _vol1, double _vol2,
    double _rho) :
    BlackScholes2D(Gen, _spot1, _spot2, _rate, _vol1, _vol2, _rho)
{
}

void BSMilstein2D::Simulate(double startTime, double endTime, size_t nbSteps)
{
    // Init Paths (2D: one path per asset)
    for (auto p : this->Paths) delete p;
    this->Paths.clear();
    this->Paths.push_back(new SinglePath(startTime, endTime, nbSteps));
    this->Paths.push_back(new SinglePath(startTime, endTime, nbSteps));

    // Simulate the 2D correlated Brownian motion using BrownianND
    // CorrelationMatrix is [[1, rho],[rho, 1]] built in BlackScholes2D
    BrownianND W(this->Generator, 2, &this->CorrelationMatrix);
    W.Simulate(startTime, endTime, nbSteps);
    SinglePath* W1_path = W.GetPath(0);
    SinglePath* W2_path = W.GetPath(1);

    double dt = (endTime - startTime) / nbSteps;

    // Initial values
    double S1 = this->Spot1;
    double S2 = this->Spot2;
    this->Paths[0]->InsertValue(S1);
    this->Paths[1]->InsertValue(S2);

    // Milstein scheme for each asset:
    // S_{i+1} = S_i + r*S_i*dt + vol*S_i*dW + 0.5*vol^2*S_i*(dW^2 - dt)
    for (size_t i = 1; i <= nbSteps; i++) {
        double t_prev = startTime + (i - 1) * dt;
        double t_curr = startTime + i * dt;

        double dW1 = W1_path->GetState(t_curr) - W1_path->GetState(t_prev);
        double dW2 = W2_path->GetState(t_curr) - W2_path->GetState(t_prev);

        S1 += this->Rate * S1 * dt
            + this->Vol1 * S1 * dW1
            + 0.5 * std::pow(this->Vol1, 2) * S1 * (std::pow(dW1, 2) - dt);

        S2 += this->Rate * S2 * dt
            + this->Vol2 * S2 * dW2
            + 0.5 * std::pow(this->Vol2, 2) * S2 * (std::pow(dW2, 2) - dt);

        this->Paths[0]->InsertValue(S1);
        this->Paths[1]->InsertValue(S2);
    }
}
