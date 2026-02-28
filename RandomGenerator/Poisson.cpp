#include "pch.h"
#include "Poisson.h"
#include <cmath>



Poisson::Poisson(double _lambda, Bernoulli _Bern) : 
    lambda(_lambda),
    Bern(_Bern)
{
}


double Poisson::Generate(PoissonAlgo algo) {

    switch (algo) {
        case PoissonAlgo::PoissonAlgo1:
            return GenerateAlgo1();
        case PoissonAlgo::PoissonAlgo2:
            return GenerateAlgo2();
        }
}


double Poisson::GenerateAlgo1() {

    double U = Bern.Generate();

    int cpt(0);
    double density_k(exp(-this->lambda));
    double cumulative(density_k);

    while (U > cumulative) {
        cpt++;
		density_k *= this->lambda / cpt;
        cumulative += density_k;
    }

    return static_cast<double>(cpt);
}


double Poisson::GenerateAlgo2() {

    return 0.0;
}



