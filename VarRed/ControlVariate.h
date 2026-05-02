#pragma once
#include "../SDE/SinglePath.h"

class ControlVariate
{
public:
    
    virtual double SimulatedValue(const std::vector<SinglePath*>& paths) const = 0;

    
    virtual double AnalyticalExpectation() const = 0;

    virtual ~ControlVariate() = default;
};

