#include "pch.h"
#include "EcuyerCombined.h"
#include <vector>


EcuyerCombined::EcuyerCombined(LinearCongruential& _firstGen, LinearCongruential& _secondGen) : PseudoGenerator(0),
    firstGen(_firstGen),
    secondGen(_secondGen)
{
};

double EcuyerCombined::Generate() {

    firstGen.Generate();
    secondGen.Generate();

    long x = static_cast<long>(firstGen.GetSeed())
        - static_cast<long>(secondGen.GetSeed());

	size_t denominator = firstGen.GetModulus() - 1;

    if (x < 0)
        x += (denominator);

    return static_cast<double>(x) / denominator;

}