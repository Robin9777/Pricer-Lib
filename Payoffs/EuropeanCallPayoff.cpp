#include "pch.h"
#include "EuropeanCallPayoff.h"
#include <cmath>

double EuropeanCallPayoff::operator()(PayOffTemplate PO_template)
{
    return std::max(PO_template.Paths.back()->GetValue() - PO_template.strike, 0.0);
}