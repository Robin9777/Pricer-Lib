#include "pch.h"
#include "RandomGenerator.h"
#include <vector>
#include <cmath>


double RandomGenerator::Mean(size_t nbSim) {

    double sum = 0.0;

    for (size_t i = 0; i < nbSim; ++i) {
        sum += this->Generate();
    }

    return sum / nbSim;
};


double RandomGenerator::Variance(size_t nbSim) {

    double mean(this->Mean(nbSim));
    double var(0);

    for (size_t i = 0; i < nbSim; i++) {
        var += std::pow(this->Generate() - mean, 2);
    }

    return var / nbSim;
}