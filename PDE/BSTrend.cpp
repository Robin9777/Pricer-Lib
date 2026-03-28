#include "pch.h"
#include "BSTrend.h"

BSTrend::BSTrend(double _rate) :
	rate(_rate)
{
}

double BSTrend::operator()(double x, double t)
{
	return rate * x;
}
