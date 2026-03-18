#include "pch.h"
#include "CppUnitTest.h"
#include <string>

#include "../RandomGenerator/Exponential.h"
#include "../RandomGenerator/LinearCongruential.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestRandomNumber
{

	TEST_CLASS(ExponentialTests)
	{
	public:

		TEST_METHOD(MeanCorrect)
		{
			size_t seed = 42, a = 1664525, c = 1013904223, m = 4294967296ULL;
			LinearCongruential unifGen(seed, a, c, m);
			
			double lambda = 2.0;
			Exponential expGen(lambda, unifGen);

			int N = 10000;
			double sum = 0.0;
			for (int i = 0; i < N; ++i)
			{
				double v = expGen.Generate();
				Assert::IsTrue(v >= 0.0, L"Exponential distribution should generate positive values.");
				sum += v;
			}
			double mean = sum / N;
			double expectedMean = 1.0 / lambda;
			// Allow some variance
			Assert::IsTrue(abs(mean - expectedMean) < 0.1, L"Mean value deviates too much from expected.");
		}

	};
}
