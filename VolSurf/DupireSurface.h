#pragma once
#include <vector>
class DupireSurface
{

private:
	std::vector<std::vector<double>> VolGrid;

public:
	double GetVariance(double spot, double time);

};

