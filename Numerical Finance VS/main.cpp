#include <iostream>
#include "../RandomGenerator/EcuyerCombined.h"


int main() {

	std::cout << "Hello World test !" << std::endl;

	LinearCongruential L1(1, 2, 3, 4);
	LinearCongruential L2(1, 2, 3, 4);
	EcuyerCombined Gen(L1, L2);

	std::cout << "Generating 10 numbers:" << std::endl;
	for (int i = 0; i < 10; ++i) {
		double v = Gen.Generate();
		std::cout << i << ": " << v << std::endl;
	}

	size_t nsim = 10000;
	std::cout << "Mean(" << nsim << ") = " << Gen.Mean(nsim) << std::endl;
	std::cout << "Variance(" << nsim << ") = " << Gen.Variance(nsim) << std::endl;

	std::cout << "Appuyez sur Entrée pour fermer..." << std::endl;
	std::cin.get();

	return 0;
}
