#include "pch.h"
#include "CppUnitTest.h"
#include "../Pricer/EuropeanMCPricer.h"
#include "../Pricer/BermudanPricer.h"
#include "../Payoffs/EuropeanCallPayoff.h"
#include "../Payoffs/EuroCallBasketPayOff.h"
#include "../SDE/BSMilstein1D.h"
#include "../SDE/BSMilstein2D.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/EcuyerCombined.h"
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
			LinearCongruential unifGenA(42, 40014, 0, 2147483563);
			LinearCongruential unifGenB(42, 40692, 0, 2147483399);
			EcuyerCombined unifGen(unifGenA, unifGenB);
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
			PriceResult res = pricer.Price();
			double price = res.price;
			
			// Log the price
			std::wstring msg = L"expected number : " + std::to_wstring(res.price) + L" +/- " + std::to_wstring(res.confidenceInterval);
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
			LinearCongruential unifGen1A(42, 40014, 0, 2147483563);
			LinearCongruential unifGen1B(42, 40692, 0, 2147483399);
			EcuyerCombined unifGen1(unifGen1A, unifGen1B);
			Normal normGen1(0.0, 1.0, unifGen1);
			BSMilstein1D bsProcess1(&normGen1, spot, rate, vol);
			EuropeanCallPayoff payoffVanilla(strike);
			EuropeanMCPricer pricerVanilla(&bsProcess1, &payoffVanilla, rate, maturity, nbSim, nbSteps);
			PriceResult resVanilla = pricerVanilla.Price();
			double priceVanilla = resVanilla.price;

			// Price with EuroCallBasketPayOff (1 asset, weight = 1.0)
			LinearCongruential unifGen2A(42, 40014, 0, 2147483563);
			LinearCongruential unifGen2B(42, 40692, 0, 2147483399);
			EcuyerCombined unifGen2(unifGen2A, unifGen2B);
			Normal normGen2(0.0, 1.0, unifGen2);
			BSMilstein1D bsProcess2(&normGen2, spot, rate, vol);
			std::vector<double> weights = { 1.0 };
			EuroCallBasketPayOff payoffBasket(strike, weights);
			EuropeanMCPricer pricerBasket(&bsProcess2, &payoffBasket, rate, maturity, nbSim, nbSteps);
			PriceResult resBasket = pricerBasket.Price();
			double priceBasket = resBasket.price;

			// Both prices should be exactly the same since seeds are identical
			std::wstring msg = L"Vanilla Price : " + std::to_wstring(resVanilla.price) + L" +/- " + std::to_wstring(resVanilla.confidenceInterval) +
			                   L" | Basket 1D Price : " + std::to_wstring(resBasket.price) + L" +/- " + std::to_wstring(resBasket.confidenceInterval);
			Logger::WriteMessage(msg.c_str());
			
			Assert::AreEqual(priceVanilla, priceBasket, 1e-2);
		}

		TEST_METHOD(TestEuroCallBasketPrice_2D)
		{
			LinearCongruential unifGenA(42, 40014, 0, 2147483563);
			LinearCongruential unifGenB(42, 40692, 0, 2147483399);
			EcuyerCombined unifGen(unifGenA, unifGenB);
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
			
			PriceResult res = pricer.Price();
			double price = res.price;
			
			std::wstring msg = L"Basket 2D Price : " + std::to_wstring(res.price) + L" +/- " + std::to_wstring(res.confidenceInterval);
			Logger::WriteMessage(msg.c_str());
			
			Assert::IsTrue(price > 0.0 && price < 15.0);
		}

		TEST_METHOD(TestBermudanCallPrice_1D_OneDate)
		{
			double spot = 100.0;
			double rate = 0.05;
			double vol = 0.2;
			double strike = 100.0;
			double maturity = 1.0;
			size_t nbSim = 2000;
			size_t nbSteps = 100;
			
			std::vector<double> exerciseDates = { maturity };

			LinearCongruential unifGen1A(42, 40014, 0, 2147483563);
			LinearCongruential unifGen1B(42, 40692, 0, 2147483399);
			EcuyerCombined unifGen1(unifGen1A, unifGen1B);
			Normal normGen1(0.0, 1.0, unifGen1);
			BSMilstein1D bsProcess1(&normGen1, spot, rate, vol);
			EuropeanCallPayoff payoffVanilla(strike);
			
			BermudanPricer pricerBermudan(&bsProcess1, &payoffVanilla, rate, maturity, nbSim, nbSteps, exerciseDates);
			PriceResult resBermudan = pricerBermudan.Price();
			double priceBermudan = resBermudan.price;

			LinearCongruential unifGen2A(42, 40014, 0, 2147483563);
			LinearCongruential unifGen2B(42, 40692, 0, 2147483399);
			EcuyerCombined unifGen2(unifGen2A, unifGen2B);
			Normal normGen2(0.0, 1.0, unifGen2);
			BSMilstein1D bsProcess2(&normGen2, spot, rate, vol);
			EuropeanMCPricer pricerEuropean(&bsProcess2, &payoffVanilla, rate, maturity, nbSim, nbSteps);
			PriceResult resEuropean = pricerEuropean.Price();
			double priceEuropean = resEuropean.price;

			std::wstring msg = L"Bermudan (1 date) : " + std::to_wstring(resBermudan.price) + L" +/- " + std::to_wstring(resBermudan.confidenceInterval) +
			                   L" | European : " + std::to_wstring(resEuropean.price) + L" +/- " + std::to_wstring(resEuropean.confidenceInterval);
			Logger::WriteMessage(msg.c_str());
			
			Assert::AreEqual(priceEuropean, priceBermudan, 1e-6);
		}

		TEST_METHOD(TestBermudanBasketPrice_2D)
		{
			LinearCongruential unifGenA(42, 40014, 0, 2147483563);
			LinearCongruential unifGenB(42, 40692, 0, 2147483399);
			EcuyerCombined unifGen(unifGenA, unifGenB);
			Normal normGen(0.0, 1.0, unifGen);
			
			double spot1 = 100.0, spot2 = 100.0;
			double rate = 0.05;
			double vol1 = 0.2, vol2 = 0.2;
			double rho = 0.5;
			BSMilstein2D bsProcess(&normGen, spot1, spot2, rate, vol1, vol2, rho);
			
			double strike = 100.0;
			std::vector<double> weights = { 0.5, 0.5 };
			EuroCallBasketPayOff payoffBasket(strike, weights);
			
			double maturity = 1.0;
			size_t nbSim = 2000;
			size_t nbSteps = 100;
			std::vector<double> exerciseDates = { 0.25, 0.5, 0.75, 1.0 };
			
			BermudanPricer pricerBermudan(&bsProcess, &payoffBasket, rate, maturity, nbSim, nbSteps, exerciseDates);
			PriceResult resBermudan = pricerBermudan.Price();
			double priceBermudan = resBermudan.price;
			
			std::wstring msg = L"Bermudan Basket 2D Price : " + std::to_wstring(resBermudan.price) + L" +/- " + std::to_wstring(resBermudan.confidenceInterval);
			Logger::WriteMessage(msg.c_str());
			
			Assert::IsTrue(priceBermudan > 0.0 && priceBermudan < 15.0);
		}
	};
}
