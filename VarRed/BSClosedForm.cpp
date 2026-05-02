#include "pch.h"
#include "BSClosedForm.h"
#include <cmath>
#include <stdexcept>

double BSClosedForm::NormalCDF(double x) const
{
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

void BSClosedForm::ComputeD1D2(double spot, double strike, double rate, double vol, double maturity, double& d1, double& d2) const
{
    const double sqrtT = std::sqrt(maturity);
    const double volSqrtT = vol * sqrtT;

    d1 = (std::log(spot / strike) + (rate + 0.5 * vol * vol) * maturity) / volSqrtT;
    d2 = d1 - volSqrtT;
}

double BSClosedForm::CallPrice(double spot, double strike, double rate, double vol, double maturity) const
{

    if (maturity <= 0.0)
        return std::max(spot - strike, 0.0);

    double d1, d2;
    ComputeD1D2(spot, strike, rate, vol, maturity, d1, d2);

    return spot * NormalCDF(d1)
        - strike * std::exp(-rate * maturity) * NormalCDF(d2);
}

double BSClosedForm::PutPrice(double spot, double strike, double rate, double vol, double maturity) const
{
    if (maturity <= 0.0)
        return std::max(strike - spot, 0.0);

    return CallPrice(spot, strike, rate, vol, maturity)
        - spot
        + strike * std::exp(-rate * maturity);
}
