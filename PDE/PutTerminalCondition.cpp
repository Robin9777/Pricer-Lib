#include "pch.h"
#include "PutTerminalCondition.h"
#include <algorithm>

PutTerminalCondition::PutTerminalCondition(double _strike) 
	: VanillaTerminalCondition(_strike)
{

}

double PutTerminalCondition::operator()(double x) {
	return std::max(this->Strike - x, 0.0)
}