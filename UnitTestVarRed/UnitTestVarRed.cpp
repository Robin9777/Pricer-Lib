#include "pch.h"
#include "CppUnitTest.h"
#include "../Pricer/EuropeanMCPricer.h"
#include "../Pricer/VarRedMCPricer.h"
#include "../Payoffs/EuroCallBasketPayOff.h"
#include "../SDE/BSMilstein2D.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/Normal.h"
#include "../VarRed/BasketGeomControlVariate.h"
#include <string>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestVarRed
{
	TEST_CLASS(UnitTestVarRed)
	{
	public:
		
		TEST_METHOD(TestVarRedBasket)
		{
			double spot1 = 100.0, spot2 = 100.0;
			double rate = 0.05;
			double vol1 = 0.2, vol2 = 0.2;
			double rho = 0.5;
			double strike = 100.0;
			double maturity = 1.0;
			size_t nbSim = 10000;
			size_t nbSteps = 100;
			
			std::vector<double> weights = { 0.5, 0.5 };
			std::vector<double> spots = { spot1, spot2 };
			std::vector<double> vols = { vol1, vol2 };
			std::vector<std::vector<double>> correlMatrix = { {1.0, rho}, {rho, 1.0} };
			
			EuroCallBasketPayOff payoffBasket(strike, weights);
			BasketGeomControlVariate cvBasket(weights, spots, correlMatrix, vols, maturity, strike, rate);
			
			// Vanilla MC
			LinearCongruential unifGen1(42, 16807, 0, 2147483647);
			Normal normGen1(0.0, 1.0, unifGen1);
			BSMilstein2D bsProcess1(&normGen1, spot1, spot2, rate, vol1, vol2, rho);
			EuropeanMCPricer pricerVanilla(&bsProcess1, &payoffBasket, rate, maturity, nbSim, nbSteps);
			PriceResult resVanilla = pricerVanilla.Price();
			
			// VarRed MC
			LinearCongruential unifGen2(42, 16807, 0, 2147483647);
			Normal normGen2(0.0, 1.0, unifGen2);
			BSMilstein2D bsProcess2(&normGen2, spot1, spot2, rate, vol1, vol2, rho);
			VarRedMCPricer pricerVarRed(&bsProcess2, &payoffBasket, &cvBasket, rate, maturity, nbSim, nbSteps);
			PriceResult resVarRed = pricerVarRed.Price();
			
			std::wstring msg = L"Vanilla CV: " + std::to_wstring(resVanilla.confidenceInterval) + L" | VarRed CV: " + std::to_wstring(resVarRed.confidenceInterval);
			Logger::WriteMessage(msg.c_str());
			
			// The control variate should significantly reduce the variance
			Assert::IsTrue(resVarRed.confidenceInterval < resVanilla.confidenceInterval);
		}
	};
}
