#include "pch.h"
#include "CppUnitTest.h"
#include <string>
#include <cmath>
#include <vector>

#include "../SDE/BSEuler1D.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/Normal.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestSDE
{
    TEST_CLASS(BSEuler1DTests)
    {
    public:

        // ----------------------------------------------------------------
        // Test 1 : S(0) = Spot  (valeur initiale dans le path)
        // ----------------------------------------------------------------
        TEST_METHOD(InitialValueIsSpot)
        {
            LinearCongruential lc(42, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            double spot = 100.0, rate = 0.05, vol = 0.2;
            BSEuler1D bs(&ng, spot, rate, vol);
            bs.Simulate(0.0, 1.0, 100);

            double s0 = bs.GetPath(0)->GetState(0.0);

            Logger::WriteMessage((std::wstring(L"S(0) généré  : ") + std::to_wstring(s0)).c_str());
            Logger::WriteMessage((std::wstring(L"S(0) attendu : ") + std::to_wstring(spot)).c_str());

            Assert::AreEqual(spot, s0, 1e-10, L"S(0) doit être égal au spot initial.");
        }

        // ----------------------------------------------------------------
        // Test 2 : Nombre de valeurs = nbSteps + 1
        // ----------------------------------------------------------------
        TEST_METHOD(PathSizeIsNbStepsPlusOne)
        {
            LinearCongruential lc(1, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            size_t nbSteps = 60;
            BSEuler1D bs(&ng, 100.0, 0.05, 0.2);
            bs.Simulate(0.0, 1.0, nbSteps);

            size_t valeur  = bs.GetPath(0)->GetAllValues().size();
            size_t attendu = nbSteps + 1;

            Logger::WriteMessage((std::wstring(L"Taille générée  : ") + std::to_wstring(valeur)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille attendue : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, valeur, L"Path size doit être nbSteps + 1.");
        }

        // ----------------------------------------------------------------
        // Test 3 : vol=0 => S(T) = S0 * exp(r*T) (déterministe)
        //          L'Euler avec vol=0 donne S_{i+1} = S_i*(1+r*dt)
        //          => S_n = S0 * (1+r*dt)^n  ≈ S0*exp(r*T) pour dt petit
        // ----------------------------------------------------------------
        TEST_METHOD(ZeroVol_PriceGrowsAtRiskFreeRate)
        {
            LinearCongruential lc(10, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            double spot = 100.0, rate = 0.05, vol = 0.0, T = 1.0;
            size_t nbSteps = 10000;  // grand N => convergence vers l'exponentielle
            BSEuler1D bs(&ng, spot, rate, vol);
            bs.Simulate(0.0, T, nbSteps);

            double sT     = bs.GetPath(0)->GetState(T);
            double dt     = T / nbSteps;
            double attendu = spot * std::pow(1.0 + rate * dt, (double)nbSteps);

            Logger::WriteMessage((std::wstring(L"S(T) généré  : ") + std::to_wstring(sT)).c_str());
            Logger::WriteMessage((std::wstring(L"S(T) attendu : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, sT, 1e-6, L"Avec vol=0, S(T) = S0*(1+r*dt)^n (schéma exact).");
        }

        // ----------------------------------------------------------------
        // Test 4 : E[S(T)] ≈ S0 * exp(r*T)  (propriété risk-neutral)
        //          (Euler est biaisé, mais convergent en loi)
        // ----------------------------------------------------------------
        TEST_METHOD(MeanPriceConvergesToForwardPrice)
        {
            size_t N = 5000;
            double spot = 100.0, rate = 0.05, vol = 0.2, T = 1.0;
            size_t nbSteps = 200;

            double sum = 0.0;
            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                BSEuler1D bs(&ng, spot, rate, vol);
                bs.Simulate(0.0, T, nbSteps);
                sum += bs.GetPath(0)->GetState(T);
            }

            double mean    = sum / N;
            double attendu = spot * std::exp(rate * T);

            Logger::WriteMessage((std::wstring(L"E[S(T)] généré  : ") + std::to_wstring(mean)).c_str());
            Logger::WriteMessage((std::wstring(L"E[S(T)] attendu : ") + std::to_wstring(attendu)).c_str());

            // Tolérance de 2 % (biais Euler + variance MC)
            Assert::IsTrue(std::abs(mean - attendu) / attendu < 0.02,
                L"E[S(T)] doit être proche de S0*exp(r*T).");
        }

        // ----------------------------------------------------------------
        // Test 5 : S(t) > 0 pour tout t  (propriété GBM)
        //          (l'Euler peut théoriquement devenir négatif, mais rarement)
        // ----------------------------------------------------------------
        TEST_METHOD(AllPathValuesPositive)
        {
            size_t N = 1000;
            double spot = 100.0, rate = 0.05, vol = 0.3, T = 1.0;
            size_t nbSteps = 252;

            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                BSEuler1D bs(&ng, spot, rate, vol);
                bs.Simulate(0.0, T, nbSteps);

                std::vector<double> vals = bs.GetPath(0)->GetAllValues();
                for (double v : vals)
                {
                    Assert::IsTrue(v > 0.0, L"Une valeur négative a été détectée dans le path.");
                }
            }
        }

        // ----------------------------------------------------------------
        // Test 6 : Ré-simulation — le path précédent est bien effacé
        // ----------------------------------------------------------------
        TEST_METHOD(ResimulationClearsPreviousPath)
        {
            LinearCongruential lc(5, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            BSEuler1D bs(&ng, 100.0, 0.05, 0.2);
            bs.Simulate(0.0, 1.0, 100);
            bs.Simulate(0.0, 1.0, 20);   // ré-simulation plus courte

            size_t valeur  = bs.GetPath(0)->GetAllValues().size();
            size_t attendu = 20 + 1;

            Logger::WriteMessage((std::wstring(L"Taille après ré-simulation : ") + std::to_wstring(valeur)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille attendue            : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, valeur, L"Après ré-simulation, le path doit avoir la nouvelle taille.");
        }
    };
}
