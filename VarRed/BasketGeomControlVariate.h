#pragma once
#include "ControlVariate.h"
class BasketGeomControlVariate :
    public ControlVariate
{
private:
	std::vector<double> Weights;
	std::vector<double> Spots;
	std::vector<std::vector<double>> correlMatrix;
	std::vector<double> Vols;

public:

	double SimulatedValue(const std::vector<SinglePath*>& paths) const override;

	double AnalyticalExpectation() const override;

};

