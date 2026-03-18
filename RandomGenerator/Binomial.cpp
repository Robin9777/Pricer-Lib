#include "pch.h"
#include "Binomial.h"


Binomial::Binomial(int _n, Bernoulli& _bern) :
	bern(_bern),
	n(_n) {

}

double Binomial::Generate() {

	unsigned int sum(0);

	for (int i = 0; i < n; i++) {
		sum += bern.Generate();
	}

	return static_cast<double>(sum);
}