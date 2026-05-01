#include "pch.h"
#include "CppUnitTest.h"
#include "../Pricer/EuropeanMCPricer.h"
#include "../Payoffs/EuropeanCallPayoff.h"
#include "../SDE/BSMilstein1D.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/Normal.h"
#include <string>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestPricer
{
	TEST_CLASS(UnitTestPricer)
	{
	public:
		
		TEST_METHOD(TestEuropeanCallPrice)
		{
			// Setup random generators
			LinearCongruential unifGen(42, 16807, 0, 2147483647);
			Normal normGen(0.0, 1.0, unifGen);
			
			// Setup Black Scholes process
			double spot = 100.0;
			double rate = 0.05;
			double vol = 0.2;
			BSMilstein1D bsProcess(&normGen, spot, rate, vol);
			
			// Setup Payoff
			double strike = 100.0;
			EuropeanCallPayoff payoff(strike);
			
			// Setup Pricer
			double maturity = 1.0;
			size_t nbSim = 10000;
			size_t nbSteps = 100;
			EuropeanMCPricer pricer(&bsProcess, &payoff, rate, maturity, nbSim, nbSteps);
			
			// Compute price
			double price = pricer.Price();
			
			// Log the price
			std::wstring msg = L"expected number : " + std::to_wstring(price);
			Logger::WriteMessage(msg.c_str());
			
			// A rough check for the price of ATM call (spot=100, K=100, r=0.05, v=0.2, T=1)
			// Expected Black-Scholes price is ~10.45
			Assert::IsTrue(price > 9.0 && price < 12.0);
		}
	};
}
