#include "pch.h"
#include "BrownianD1.h"
#include "RandomProcess.h"
#include "../SDE/SinglePath.h"
#include <cmath>

Brownian1D::Brownian1D(RandomGenerator* Gen) : RandomProcess(Gen, 1)
{
}

void Brownian1D::Simulate(double startTime, double endTime, size_t nbSteps)
{
    // Clean up existing paths
    for (auto p : this->Paths) {
        delete p;
    }
    this->Paths.clear();

    // Create new path
    this->Paths.push_back(new SinglePath(startTime, endTime, nbSteps));

    double dt = (endTime - startTime) / nbSteps;
    double sqrt_dt = std::sqrt(dt);

    double current_value = 0.0;
    this->Paths[0]->InsertValue(current_value);

    for (size_t i = 1; i <= nbSteps; i++) {
        double sim = this->Generator->Generate();
        double increment = sim * sqrt_dt;
        current_value += increment;
        this->Paths[0]->InsertValue(current_value);
    }
}
