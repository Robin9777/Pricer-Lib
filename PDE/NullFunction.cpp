#include "pch.h"
#include "NullFunction.h"

NullFunction::NullFunction()
{
}

double NullFunction::operator()(double x, double t) {
	return 0.0;
}