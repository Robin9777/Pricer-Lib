#include "pch.h"
#include "BSMilstein1D.h"
#include "BrownianD1.h"
#include "../SDE/SinglePath.h"
#include <cmath>

BSMilstein1D::BSMilstein1D(RandomGenerator* Gen, double _spot, double _rate, double _vol) :
	BlackScholes1D(Gen, _spot, _rate, _vol)
{
}

void BSMilstein1D::Simulate(double startTime, double endTime, size_t nbSteps)
{
	// Init Path
	for (auto p : this->Paths) delete p;
	this->Paths.clear();
	this->Paths.push_back(new SinglePath(startTime, endTime, nbSteps));

	// Simulate the brownian motion
	Brownian1D W(this->Generator);
	W.Simulate(startTime, endTime, nbSteps);
	SinglePath* W_path = W.GetPath(0);

	double dt = (endTime - startTime) / nbSteps;
	double currentvalue = this->Spot;
	this->Paths[0]->InsertValue(currentvalue);

	// Looping
	for (size_t i = 1; i <= nbSteps; i++) {
		double t_prev = startTime + (i - 1) * dt;
		double t_curr = startTime + i * dt;
		double dW = W_path->GetState(t_curr) - W_path->GetState(t_prev);

		// currentvalue += (this->Rate - 1/2 * std::pow(this->Vol, 2)) * currentvalue * dt
		// 	+ 1/2 * std::pow(this->Vol, 2) * currentvalue * std::pow(dW, 2);
		// Milstein : S_{i+1} = S_i + r*S_i*dt + vol*S_i*dW + 0.5*vol^2*S_i*(dW^2 - dt)
		currentvalue += this->Rate * currentvalue * dt
			+ this->Vol * currentvalue * dW
			+ 0.5 * std::pow(this->Vol, 2) * currentvalue * (std::pow(dW, 2) - dt);

		this->Paths[0]->InsertValue(currentvalue);
	}
}
