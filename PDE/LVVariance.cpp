#include "pch.h"
#include "LVVariance.h"

double LVVariance::operator()(double x, double t)
{

	double variance = this->Surface->GetVariance(x, t);
    return 0.5 * variance * x*x;
}
