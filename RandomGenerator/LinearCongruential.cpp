#include "pch.h"
#include "LinearCongruential.h"


LinearCongruential::LinearCongruential(size_t _seed, size_t _multiplier, size_t _increment, size_t _modulus) :
	PseudoGenerator(_seed),
	multiplier(_multiplier),
	increment(_increment),
	modulus(_modulus)
{
	
}


double LinearCongruential::Generate() {
	this->SetSeed((this->multiplier * this->GetSeed() + this->increment) % this->modulus);
	return static_cast<double>(this->GetSeed()) / this->modulus;
}

size_t LinearCongruential::GetModulus() const {
	return this->modulus;
}
