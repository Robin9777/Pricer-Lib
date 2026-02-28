#pragma once
#include "PseudoGenerator.h"
#include "LinearCongruential.h"


class EcuyerCombined :
    public PseudoGenerator
{

protected:
    LinearCongruential& firstGen;
    LinearCongruential& secondGen;

public:
    EcuyerCombined() = default;
    EcuyerCombined(LinearCongruential& _firstGen, LinearCongruential& _secondGen);

    virtual double Generate();

};

