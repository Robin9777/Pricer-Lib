#include "pch.h"
#include "RandomProcess.h"

RandomProcess::RandomProcess(RandomGenerator* Gen, int dim) :
	Generator(Gen),
	Dimension(dim)
{
}

SinglePath* RandomProcess::GetPath(int dimension)
{
	return this->Paths[dimension];
}
