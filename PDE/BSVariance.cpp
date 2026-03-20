#include "pch.h"
#include "BSVariance.h"

BSVariance::BSVariance(double _sigma) :
	Sigma(_sigma)
{
}

double BSVariance::operator()(double x, double t)
{
	return Sigma*Sigma*x*x;
}
