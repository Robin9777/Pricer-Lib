#pragma once
#include <vector>
#include "../SDE/SinglePath.h"
class PayOff
{

public:
	PayOff() = default;
	virtual double operator()(const std::vector<SinglePath*>& Paths) const = 0;
	virtual ~PayOff() = default;
};

