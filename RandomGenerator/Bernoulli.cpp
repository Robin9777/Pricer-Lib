#include "pch.h"
#include "Bernoulli.h"
#include "UniformGenerator.h"


Bernoulli::Bernoulli(const double& _p, UniformGenerator& _Ugen) : Ugen(_Ugen), p(_p) {

}

double Bernoulli::Generate() {
	return (Ugen.Generate() < p) ? 1 : 0;
}