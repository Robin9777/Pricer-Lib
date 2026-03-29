#include "pch.h"
#include "BSEuler1D.h"
#include "BrownianD1.h"

BSEuler1D::BSEuler1D(RandomGenerator* Gen, double _spot, double _rate, double _vol) :
	BlackScholes1D(Gen, _spot, _rate, _vol)
{
}

void BSEuler1D::Simulate(double startTime, double endTime, size_t nbSteps)
{
	
	// Simulate the brownian motion
	Brownian1D W(this->Generator);
	W.Simulate(startTime, endTime, nbSteps);
	SinglePath* W_path = W.GetPath(0);

	// Initialize
	double dt = (endTime - startTime) / nbSteps;
	this->Paths[0]->InsertValue(this->Spot);
	double currentvalue(this->Paths[0]->GetState(startTime));

	// looping 
	for (size_t i = 1; i <= nbSteps; i++) {

		currentvalue += this->Rate * dt + this->Vol * (W_path->GetState(i) - W_path->GetState(i - 1));

		this->Paths[0]->InsertValue(currentvalue);

	}

}