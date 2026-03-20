#include "pch.h"
#include "CallTopBoundary.h"

CallTopBoundary::CallTopBoundary(double _Smax, double _Strike) : 
	Smax(_Smax),
	Strike(_Strike)
{
}

double CallTopBoundary::operator()(double t)
{
	return this->Smax - this->Strike;
}
