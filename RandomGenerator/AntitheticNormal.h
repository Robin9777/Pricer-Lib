#pragma once
#include "RandomGenerator.h"
#include "Normal.h"
#include <vector>

class AntitheticNormal : public RandomGenerator
{
    Normal&              inner;
    std::vector<double>  buffer;
    size_t               idx;
    bool                 antithetic;

public:
    explicit AntitheticNormal(Normal& n);

    double Generate() override;

    void SetAntithetic(bool flag);
    void ResetBuffer();
    void ResetIndex();
};
