#include "pch.h"
#include "CallBottomBoundary.h"

CallBottomBoundary::CallBottomBoundary(double _Smin, double _Strike) :
	Smin(_Smin),
	Strike(_Strike)
{
}

double CallBottomBoundary::operator()(double t)
{
	return 0.0;
}
