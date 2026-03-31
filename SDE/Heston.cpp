#include "pch.h"
#include "Heston.h"
#include "BrownianND.h"
#include "../SDE/SinglePath.h"
#include <cmath>

Heston::Heston(RandomGenerator* Gen, double _spot, double _initVariance, double _mu, double _theta, double _kappa, double _sigma, double _rho) :
	RandomProcess(Gen, 2),
	Spot(_spot),
	InitVariance(_initVariance),
	Mu(_mu),
	Theta(_theta),
	Kappa(_kappa),
	Sigma(_sigma),
	Rho(_rho)
{
	CorrelationMatrix = { {1.0, _rho}, {_rho, 1.0} };
}

void Heston::Simulate(double startTime, double endTime, size_t nbSteps)
{

	// Init Path
	for (auto p : this->Paths) delete p;
	this->Paths.clear();
	this->Paths.push_back(new SinglePath(startTime, endTime, nbSteps));

	BrownianND W(this->Generator, 2, &this->CorrelationMatrix);
	W.Simulate(startTime, endTime, nbSteps);
	SinglePath* W1_path = W.GetPath(0);
	SinglePath* W2_path = W.GetPath(1);

	double dt = (endTime - startTime) / nbSteps;

	// Intiate process
	double S = this->Spot;
	double V = this->InitVariance;
	this->Paths[0]->InsertValue(S);
	this->Paths[1]->InsertValue(V);

	// Heston Euler process
	for (size_t i = 1; i <= nbSteps; i++) {

		double t_prev = startTime + (i - 1) * dt;
		double t_curr = startTime + i * dt;

		double dW1 = W1_path->GetState(t_curr) - W1_path->GetState(t_prev);
		double dW2 = W2_path->GetState(t_curr) - W2_path->GetState(t_prev);

		V += this->Kappa * (this->Theta - V) * dt + this->Sigma * std::sqrt(V) * dW2;
		S += this->Mu * S * dt + std::sqrt(V) * S * dW1;

		this->Paths[0]->InsertValue(S);
		this->Paths[1]->InsertValue(V);

	}
}
