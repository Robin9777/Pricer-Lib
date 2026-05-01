#pragma once
#include <vector>
class SinglePath
{

protected:
	std::vector<double> Values;
	double StartTime;
	double EndTime;
	size_t NbSteps;

public:
	SinglePath(double _start, double _end, size_t _nbSteps);
	void InsertValue(double val);

	// Getter
	double GetState(double time);
	const std::vector<double>& GetAllValues() const;

};

