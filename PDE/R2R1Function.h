#pragma once
class R2R1Function
{

public:
	R2R1Function() = default;
	virtual double operator()(double x, double t) = 0;
};