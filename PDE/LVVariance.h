#pragma once
#include "R2R1Function.h"
#include "../VolSurf/DupireSurface.h"

class LVVariance :
    public R2R1Function
{

private:
	DupireSurface* Surface;

public:
	double operator()(double x, double t) override;
};

