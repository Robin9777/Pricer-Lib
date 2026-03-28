#include "pch.h"
#include "SinglePath.h"

SinglePath::SinglePath(double _start, double _end, size_t _nbSteps) :
	StartTime(_start),
	EndTime(_end),
	NbSteps(_nbSteps)
{
}

void SinglePath::InsertValue(double val)
{
	Values.push_back(val);
}

double SinglePath::GetState(double time)
{

    double ratio = (time - StartTime) / (EndTime - StartTime);
    size_t index = static_cast<size_t>(ratio * NbSteps);

    if (index >= Values.size())
        index = Values.size() - 1;

    return Values[index];
}

std::vector<double> SinglePath::GetAllValues()
{
	return this->Values;
}
