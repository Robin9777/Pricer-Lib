#pragma once
#include "PayOffTemplate.h"
class PayOff
{

public:
	PayOff() = default;
	virtual double operator()(PayOffTemplate PO_template)=0;
};

