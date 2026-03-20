#include "pch.h"
#include "BSActualization.h"
#include <cmath>

BSActualization::BSActualization(double _Rate) :
	Rate(_Rate)
{
}

double BSActualization::operator()(double x, double t)
{
	return std::exp(-x*t);
}
