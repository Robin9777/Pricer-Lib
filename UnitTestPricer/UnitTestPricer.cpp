#include "pch.h"
#include "CppUnitTest.h"
#include "../Pricer/EuropeanMCPricer.h"
#include "../Payoffs/EuropeanCallPayoff.h"
#include "../Payoffs/EuroCallBasketPayOff.h"
#include "../SDE/BSMilstein1D.h"
#include "../SDE/BSMilstein2D.h"
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

		TEST_METHOD(TestEuroCallBasketPrice_1D_Equivalent)
		{
			double spot = 100.0;
			double rate = 0.05;
			double vol = 0.2;
			double strike = 100.0;
			double maturity = 1.0;
			size_t nbSim = 10000;
			size_t nbSteps = 100;

			// Price with standard EuropeanCallPayoff
			LinearCongruential unifGen1(42, 16807, 0, 2147483647);
			Normal normGen1(0.0, 1.0, unifGen1);
			BSMilstein1D bsProcess1(&normGen1, spot, rate, vol);
			EuropeanCallPayoff payoffVanilla(strike);
			EuropeanMCPricer pricerVanilla(&bsProcess1, &payoffVanilla, rate, maturity, nbSim, nbSteps);
			double priceVanilla = pricerVanilla.Price();

			// Price with EuroCallBasketPayOff (1 asset, weight = 1.0)
			LinearCongruential unifGen2(42, 16807, 0, 2147483647);
			Normal normGen2(0.0, 1.0, unifGen2);
			BSMilstein1D bsProcess2(&normGen2, spot, rate, vol);
			std::vector<double> weights = { 1.0 };
			EuroCallBasketPayOff payoffBasket(strike, weights);
			EuropeanMCPricer pricerBasket(&bsProcess2, &payoffBasket, rate, maturity, nbSim, nbSteps);
			double priceBasket = pricerBasket.Price();

			// Both prices should be exactly the same since seeds are identical
			std::wstring msg = L"Vanilla Price : " + std::to_wstring(priceVanilla) + L" | Basket 1D Price : " + std::to_wstring(priceBasket);
			Logger::WriteMessage(msg.c_str());
			
			Assert::AreEqual(priceVanilla, priceBasket, 1e-6);
		}

		TEST_METHOD(TestEuroCallBasketPrice_2D)
		{
			LinearCongruential unifGen(42, 16807, 0, 2147483647);
			Normal normGen(0.0, 1.0, unifGen);
			
			double spot1 = 100.0;
			double spot2 = 100.0;
			double rate = 0.05;
			double vol1 = 0.2;
			double vol2 = 0.2;
			double rho = 0.5;
			
			BSMilstein2D bsProcess2D(&normGen, spot1, spot2, rate, vol1, vol2, rho);
			
			double strike = 100.0;
			std::vector<double> weights = { 0.5, 0.5 };
			EuroCallBasketPayOff payoffBasket(strike, weights);
			
			double maturity = 1.0;
			size_t nbSim = 10000;
			size_t nbSteps = 100;
			EuropeanMCPricer pricer(&bsProcess2D, &payoffBasket, rate, maturity, nbSim, nbSteps);
			
			double price = pricer.Price();
			
			std::wstring msg = L"Basket 2D Price : " + std::to_wstring(price);
			Logger::WriteMessage(msg.c_str());
			
			Assert::IsTrue(price > 0.0 && price < 15.0);
		}
	};
}
