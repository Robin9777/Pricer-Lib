#include <iostream>
#include "LinearCongruential.h"


int main() {

	std::cout << "Hello World test !\n";

	LinearCongruential Gen = LinearCongruential(
		27,
		17,
		43,
		100
	);

	std::cout << "test passed.\n";

	// Petit test : générer quelques nombres, puis calculer moyenne/variance

	std::cout << "Generating 10 numbers:\n";
	for (int i = 0; i < 10; ++i) {
		double v = Gen.Generate();
		std::cout << i << ": " << v << std::endl;
	}

	size_t nsim = 10000;
	std::cout << "Mean(" << nsim << ") = " << Gen.Mean(nsim) << std::endl;
	std::cout << "Variance(" << nsim << ") = " << Gen.Variance(nsim) << std::endl;

	return 0;

}
