#include <iostream>
#include "../RandomGenerator/LinearCongruential.h"


int main() {

	std::cout << "Hello World test !" << std::endl;

	LinearCongruential Gen(
		27u,
		17u,
		43u,
		100u
	);

	std::cout << "test passed." << std::endl;

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
