#include "pch.h"
#include "HeadTail.h"
#include "UniformGenerator.h"


HeadTail::HeadTail(UniformGenerator& _Ugen) : Ugen(_Ugen) {

}

double HeadTail::Generate() {

	double uni_number;
	uni_number = Ugen.Generate();

	return (uni_number < 0.5) ? 0 : 1;

}