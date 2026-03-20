#include "pch.h"
#include "PutBottomBoundary.h"

PutBottomBoundary::PutBottomBoundary(double _Smin, double _Strike) :
	Smin(_Smin),
	Strike(_Strike)
{
}

double PutBottomBoundary::operator()(double t)
{
	return this->Strike;
}
