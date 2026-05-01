#include "pch.h"
#include "EuropeanCallPayoff.h"
#include <cmath>

EuropeanCallPayoff::EuropeanCallPayoff(double _strike) : strike(_strike)
{
}

double EuropeanCallPayoff::operator()(const std::vector<SinglePath*>& Paths) const
{
    
    return std::max(Paths[0]->GetAllValues().back() - strike, 0.0);
}