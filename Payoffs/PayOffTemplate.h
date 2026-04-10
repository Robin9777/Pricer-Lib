#pragma once
#include "../SDE/SinglePath.h"
#include <vector>

struct PayOffTemplate
{
	
	std::vector<SinglePath*> Paths;
	double strike;
	double maturity;

};