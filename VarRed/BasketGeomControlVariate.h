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
	double maturity;
	double strike;
	double rate;

public:
	BasketGeomControlVariate(
		const std::vector<double>& _weights,
		const std::vector<double>& _spots,
		const std::vector<std::vector<double>>& _correlMatrix,
		const std::vector<double>& _vols,
		double _maturity,
		double _strike,
		double _rate
	);

	double SimulatedValue(const std::vector<SinglePath*>& paths) const override;

	double AnalyticalExpectation() const override;

};

