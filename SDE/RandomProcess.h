#pragma once
#include "../RandomGenerator/RandomGenerator.h"
#include <vector>
#include "../SDE/SinglePath.h"
class RandomProcess
{

protected:
	RandomGenerator* Generator;
	std::vector<SinglePath*> Paths;
	int Dimension;

public:
	RandomProcess(RandomGenerator* Gen, int dim);
	virtual void Simulate(double startTime, double endTime, size_t nbSteps) = 0;
	SinglePath* GetPath(int dimension = 0);

};

