#include "pch.h"
#include "PutTopBoundary.h"

PutTopBoundary::PutTopBoundary(double _Smax, double _Strike) :
	Smax(_Smax),
	Strike(_Strike)
{
}

double PutTopBoundary::operator()(double t)
{
	return 0;
}
