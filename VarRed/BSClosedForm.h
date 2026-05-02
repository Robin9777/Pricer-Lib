#pragma once
class BSClosedForm
{
private:
	double NormalCDF(double x) const;
	void ComputeD1D2(double spot, double strike, double rate, double vol, double maturity, double& d1, double& d2) const;


public:
	BSClosedForm() = default;
	double CallPrice(double spot, double strike, double rate, double vol, double maturity) const;
	double PutPrice(double spot, double strike, double rate, double vol, double maturity) const;
};

