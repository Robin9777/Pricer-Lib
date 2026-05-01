#include "pch.h"
#include "EuroCallBasketPayOff.h"
#include <numeric>
#include <cmath>

EuroCallBasketPayOff::EuroCallBasketPayOff(double _strike, const std::vector<double>& _weights)
    : strike(_strike), weights(_weights)
{
}

double EuroCallBasketPayOff::operator()(const std::vector<SinglePath*>& Paths) const
{
	double basketValue = 0.0;
	for (size_t i = 0; i < weights.size(); ++i) {
		basketValue += weights[i] * Paths[i]->GetAllValues().back();
	}
	return std::max(0.0, basketValue - strike);
}
