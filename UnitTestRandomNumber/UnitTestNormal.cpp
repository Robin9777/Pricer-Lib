#include "pch.h"
#include "CppUnitTest.h"
#include <string>
#include <cmath>

#include "../RandomGenerator/Normal.h"
#include "../RandomGenerator/LinearCongruential.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestRandomNumber
{

	TEST_CLASS(NormalTests)
	{
	public:

		TEST_METHOD(MeanAndVarianceCorrect)
		{
			size_t seed = 42, a = 1664525, c = 1013904223, m = 4294967296ULL;
			LinearCongruential unifGen(seed, a, c, m);
			
			double mu = 5.0;
			double sigma = 2.0;
			Normal normalGen(mu, sigma, unifGen);

			int N = 10000;
			double sum = 0.0;
			double sumSq = 0.0;
			for (int i = 0; i < N; ++i)
			{
				double v = normalGen.Generate();
				sum += v;
				sumSq += v * v;
			}
			double mean = sum / N;
			double expectedMean = mu;
			double var = (sumSq / N) - (mean * mean);
			double expectedVar = sigma * sigma;

			// Allow some variance
			Assert::IsTrue(std::abs(mean - expectedMean) < 0.1, L"Mean value deviates too much from expected.");
			Assert::IsTrue(std::abs(var - expectedVar) < 0.2, L"Variance deviates too much from expected.");
		}

	};
}
