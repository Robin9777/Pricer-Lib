#include "pch.h"
#include "AntitheticNormal.h"

AntitheticNormal::AntitheticNormal(Normal& n)
    : inner(n), idx(0), antithetic(false)
{
}

double AntitheticNormal::Generate()
{
    if (!antithetic) {
        double z = inner.Generate();
        buffer.push_back(z);
        return z;
    }
    return -buffer[idx++];
}

void AntitheticNormal::SetAntithetic(bool flag) { antithetic = flag; }
void AntitheticNormal::ResetBuffer()             { buffer.clear(); idx = 0; }
void AntitheticNormal::ResetIndex()              { idx = 0; }
