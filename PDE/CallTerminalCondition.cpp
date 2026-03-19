#include "pch.h"
#include "CallTerminalCondition.h"
#include <algorithm>

CallTerminalCondition::CallTerminalCondition(double _strike)
	: VanillaTerminalCondition(_strike)
{
}

double CallTerminalCondition::operator()(double x) {
	return std::max(x - this->Strike, 0.0);
}