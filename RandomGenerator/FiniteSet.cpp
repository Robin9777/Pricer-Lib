#include "pch.h"
#include "FiniteSet.h"


FiniteSet::FiniteSet(std::vector<double> _values, Bernoulli _Bern) :
	values(_values),
	Bern(_Bern)

{
}

double FiniteSet::Generate() {

    double sum = 0.0;
    for (double v : values) {
        sum += v;
    }

    double U = Bern.Generate();

    double cumulative = 0.0;

    for (size_t k = 0; k < values.size(); ++k) {
        cumulative += values[k] / sum;
        if (U < cumulative) {
            return static_cast<double>(k);
        }
    }

    return static_cast<double>(values.size() - 1);
}