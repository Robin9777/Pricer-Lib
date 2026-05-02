#include "pch.h"
#include "BasketGeomControlVariate.h"
#include <cmath>
#include "BSClosedForm.h"

BasketGeomControlVariate::BasketGeomControlVariate(
	const std::vector<double>& _weights,
	const std::vector<double>& _spots,
	const std::vector<std::vector<double>>& _correlMatrix,
	const std::vector<double>& _vols,
	double _maturity,
	double _strike,
	double _rate
) : Weights(_weights), Spots(_spots), correlMatrix(_correlMatrix), Vols(_vols), maturity(_maturity), strike(_strike), rate(_rate) {}

double BasketGeomControlVariate::SimulatedValue(const std::vector<SinglePath*>& paths) const
{
	double logGeomSum = 0.0;
	// for 
	for (size_t i = 0; i < Weights.size(); ++i) {
		double S_i_T = paths[i]->GetAllValues().back();
		logGeomSum += Weights[i] * std::log(S_i_T);
	}
	double geomBasket = std::exp(logGeomSum);

	return std::max(geomBasket - strike, 0.0);
}

double BasketGeomControlVariate::AnalyticalExpectation() const
{

	// Init BlackScholes Option
	BSClosedForm bs;

	std::vector<std::vector<double>> covMatrix(Spots.size(), std::vector<double>(Spots.size(), 0.0));
	// Compute the covariance matrix
	for (size_t i = 0; i < Spots.size(); ++i) {
		for (size_t j = 0; j < Spots.size(); ++j) {
			covMatrix[i][j] = Vols[i] * Vols[j] * correlMatrix[i][j];
		}
	}

	// Compute the variance of the geometric average
	double varGeom = 0.0;
	for (size_t i = 0; i < Spots.size(); ++i) {
		for (size_t j = 0; j < Spots.size(); ++j) {
			varGeom += Weights[i] * Weights[j] * covMatrix[i][j];
		}
	}

	// Geometric spot price
	double geomSpot = 1.0;
	for (size_t i = 0; i < Spots.size(); ++i) {
		geomSpot *= std::pow(Spots[i], Weights[i]);
	}


    return bs.CallPrice(geomSpot, strike, rate, varGeom, maturity);
}
