#include "pch.h"
#include "CppUnitTest.h"
#include <string>

#include "../RandomGenerator/EcuyerCombined.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestRandomNumber
{

	TEST_CLASS(EcuyerCombinedTests)
	{
	public:

		TEST_METHOD(ValueIn01)
		{
			size_t seed_1 = 42, a_1 = 1664525, c_1 = 1013904223, m_1 = 4294967296ULL;
			size_t seed_2 = 42, a_2 = 1664525, c_2 = 1013904223, m_2 = 4294967296ULL;

			LinearCongruential l1(seed_1, a_1, c_1, m_1);
			LinearCongruential l2(seed_2, a_2, c_2, m_2);

			EcuyerCombined gen(l1, l2);

			int N = 1000;
			for (int i = 0; i < N; ++i)
			{
				double v = gen.Generate();
				Assert::IsTrue(v >= 0.0 && v < 1.0,
					L"Une valeur générée est en dehors de l'intervalle [0, 1).");
			}
		}



	};
}