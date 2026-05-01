#include "pch.h"
#include "CppUnitTest.h"
#include "../Payoffs/EuropeanCallPayoff.h"
#include "../SDE/SinglePath.h"
#include <vector>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestPayoffs
{
	TEST_CLASS(UnitTestPayoffs)
	{
	public:
		
		TEST_METHOD(TestEuropeanCallOTM)
		{
			EuropeanCallPayoff payoff(100.0);
			SinglePath path(0.0, 1.0, 1);
			path.InsertValue(100.0);
			path.InsertValue(90.0); // S_T = 90 < 100
			std::vector<SinglePath*> paths = { &path };
			Assert::AreEqual(0.0, payoff(paths));
		}

		TEST_METHOD(TestEuropeanCallITM)
		{
			EuropeanCallPayoff payoff(100.0);
			SinglePath path(0.0, 1.0, 1);
			path.InsertValue(100.0);
			path.InsertValue(110.0); // S_T = 110 > 100
			std::vector<SinglePath*> paths = { &path };
			Assert::AreEqual(10.0, payoff(paths));
		}
		
		TEST_METHOD(TestEuropeanCallATM)
		{
			EuropeanCallPayoff payoff(100.0);
			SinglePath path(0.0, 1.0, 1);
			path.InsertValue(100.0);
			path.InsertValue(100.0); // S_T = 100
			std::vector<SinglePath*> paths = { &path };
			Assert::AreEqual(0.0, payoff(paths));
		}
	};
}
