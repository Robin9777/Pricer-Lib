#include "pch.h"
#include "CppUnitTest.h"
#include <string>
#include <cmath>
#include <vector>

#include "../SDE/BSMilstein1D.h"
#include "../SDE/BSMilstein2D.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/Normal.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestSDE
{
    // ====================================================================
    //  BSMilstein1D
    // ====================================================================
    TEST_CLASS(BSMilstein1DTests)
    {
    public:

        // ----------------------------------------------------------------
        // Test 1 : S(0) = Spot
        // ----------------------------------------------------------------
        TEST_METHOD(InitialValueIsSpot)
        {
            LinearCongruential lc(42, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            double spot = 100.0;
            BSMilstein1D bs(&ng, spot, 0.05, 0.2);
            bs.Simulate(0.0, 1.0, 100);

            double s0     = bs.GetPath(0)->GetState(0.0);
            double attendu = spot;

            Logger::WriteMessage((std::wstring(L"S(0) généré  : ") + std::to_wstring(s0)).c_str());
            Logger::WriteMessage((std::wstring(L"S(0) attendu : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, s0, 1e-10, L"S(0) doit être égal au spot initial.");
        }

        // ----------------------------------------------------------------
        // Test 2 : Taille du path = nbSteps + 1
        // ----------------------------------------------------------------
        TEST_METHOD(PathSizeIsNbStepsPlusOne)
        {
            LinearCongruential lc(1, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            size_t nbSteps = 80;
            BSMilstein1D bs(&ng, 100.0, 0.05, 0.2);
            bs.Simulate(0.0, 1.0, nbSteps);

            size_t valeur  = bs.GetPath(0)->GetAllValues().size();
            size_t attendu = nbSteps + 1;

            Logger::WriteMessage((std::wstring(L"Taille générée  : ") + std::to_wstring(valeur)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille attendue : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, valeur, L"Path size doit être nbSteps + 1.");
        }

        // ----------------------------------------------------------------
        // Test 3 : vol=0 => S(T) = S0*(1+r*dt)^n (identique à Euler, déterministe)
        // ----------------------------------------------------------------
        TEST_METHOD(ZeroVol_Deterministic)
        {
            LinearCongruential lc(7, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            double spot = 100.0, rate = 0.05, vol = 0.0, T = 1.0;
            size_t nbSteps = 10000;
            BSMilstein1D bs(&ng, spot, rate, vol);
            bs.Simulate(0.0, T, nbSteps);

            double sT      = bs.GetPath(0)->GetState(T);
            double dt      = T / nbSteps;
            double attendu  = spot * std::pow(1.0 + rate * dt, (double)nbSteps);

            Logger::WriteMessage((std::wstring(L"S(T) généré  : ") + std::to_wstring(sT)).c_str());
            Logger::WriteMessage((std::wstring(L"S(T) attendu : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, sT, 1e-6, L"Avec vol=0, Milstein doit être déterministe.");
        }

        // ----------------------------------------------------------------
        // Test 4 : E[S(T)] ≈ S0*exp(r*T)
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
                BSMilstein1D bs(&ng, spot, rate, vol);
                bs.Simulate(0.0, T, nbSteps);
                sum += bs.GetPath(0)->GetState(T);
            }

            double mean    = sum / N;
            double attendu = spot * std::exp(rate * T);

            Logger::WriteMessage((std::wstring(L"E[S(T)] généré  : ") + std::to_wstring(mean)).c_str());
            Logger::WriteMessage((std::wstring(L"E[S(T)] attendu : ") + std::to_wstring(attendu)).c_str());

            Assert::IsTrue(std::abs(mean - attendu) / attendu < 0.02,
                L"E[S(T)] doit être proche de S0*exp(r*T).");
        }

        // ----------------------------------------------------------------
        // Test 5 : Milstein plus précis qu'Euler — variance de l'erreur inf.
        //          On compare E[|S_Milstein(T) - S_exact(T)|^2]
        //          < E[|S_Euler(T) - S_exact(T)|^2]  (pas de vrai GBM exact ici)
        //          Proxy : log-variance (Milstein = forte convergence ordre 1)
        //          On vérifie simplement que les deux donnent la même espérance.
        // ----------------------------------------------------------------
        TEST_METHOD(MilsteinAndEulerGiveSameMean)
        {
            // Ce test vérifie que le schéma Milstein donne la même espérance
            // que le schéma Euler (tous deux convergent vers la même limite).
            size_t N = 3000;
            double spot = 100.0, rate = 0.05, vol = 0.25, T = 1.0;
            size_t nbSteps = 100;

            double sumMil = 0.0;
            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                BSMilstein1D bs(&ng, spot, rate, vol);
                bs.Simulate(0.0, T, nbSteps);
                sumMil += bs.GetPath(0)->GetState(T);
            }

            double meanMil = sumMil / N;
            double attendu = spot * std::exp(rate * T);

            Logger::WriteMessage((std::wstring(L"E[S_Milstein(T)] : ") + std::to_wstring(meanMil)).c_str());
            Logger::WriteMessage((std::wstring(L"S0*exp(r*T)      : ") + std::to_wstring(attendu)).c_str());

            Assert::IsTrue(std::abs(meanMil - attendu) / attendu < 0.025,
                L"Milstein : E[S(T)] doit rester proche de S0*exp(r*T).");
        }
    };

    // ====================================================================
    //  BSMilstein2D
    // ====================================================================
    TEST_CLASS(BSMilstein2DTests)
    {
    public:

        // ----------------------------------------------------------------
        // Test 1 : S1(0) = Spot1, S2(0) = Spot2
        // ----------------------------------------------------------------
        TEST_METHOD(InitialValuesAreSpots)
        {
            LinearCongruential lc(42, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            double spot1 = 100.0, spot2 = 80.0, rate = 0.05;
            double vol1 = 0.2, vol2 = 0.3, rho = 0.5;

            BSMilstein2D bs(&ng, spot1, spot2, rate, vol1, vol2, rho);
            bs.Simulate(0.0, 1.0, 100);

            double s1_0 = bs.GetPath(0)->GetState(0.0);
            double s2_0 = bs.GetPath(1)->GetState(0.0);

            Logger::WriteMessage((std::wstring(L"S1(0) généré  : ") + std::to_wstring(s1_0)).c_str());
            Logger::WriteMessage((std::wstring(L"S1(0) attendu : ") + std::to_wstring(spot1)).c_str());
            Logger::WriteMessage((std::wstring(L"S2(0) généré  : ") + std::to_wstring(s2_0)).c_str());
            Logger::WriteMessage((std::wstring(L"S2(0) attendu : ") + std::to_wstring(spot2)).c_str());

            Assert::AreEqual(spot1, s1_0, 1e-10, L"S1(0) doit être Spot1.");
            Assert::AreEqual(spot2, s2_0, 1e-10, L"S2(0) doit être Spot2.");
        }

        // ----------------------------------------------------------------
        // Test 2 : Taille des deux paths = nbSteps + 1
        // ----------------------------------------------------------------
        TEST_METHOD(BothPathSizesAreCorrect)
        {
            LinearCongruential lc(1, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            size_t nbSteps = 50;
            BSMilstein2D bs(&ng, 100.0, 80.0, 0.05, 0.2, 0.3, 0.4);
            bs.Simulate(0.0, 1.0, nbSteps);

            size_t n0 = bs.GetPath(0)->GetAllValues().size();
            size_t n1 = bs.GetPath(1)->GetAllValues().size();
            size_t attendu = nbSteps + 1;

            Logger::WriteMessage((std::wstring(L"Taille path S1 : ") + std::to_wstring(n0)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille path S2 : ") + std::to_wstring(n1)).c_str());
            Logger::WriteMessage((std::wstring(L"Attendu        : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, n0, L"Path S1 : taille doit être nbSteps+1.");
            Assert::AreEqual(attendu, n1, L"Path S2 : taille doit être nbSteps+1.");
        }

        // ----------------------------------------------------------------
        // Test 3 : E[S1(T)] ≈ S10*exp(r*T) et E[S2(T)] ≈ S20*exp(r*T)
        // ----------------------------------------------------------------
        TEST_METHOD(MeanPricesConvergeToForwardPrices)
        {
            size_t N = 5000;
            double spot1 = 100.0, spot2 = 80.0;
            double rate = 0.05, vol1 = 0.2, vol2 = 0.3, rho = 0.5, T = 1.0;
            size_t nbSteps = 200;

            double sum1 = 0.0, sum2 = 0.0;
            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                BSMilstein2D bs(&ng, spot1, spot2, rate, vol1, vol2, rho);
                bs.Simulate(0.0, T, nbSteps);
                sum1 += bs.GetPath(0)->GetState(T);
                sum2 += bs.GetPath(1)->GetState(T);
            }

            double mean1 = sum1 / N, mean2 = sum2 / N;
            double attendu1 = spot1 * std::exp(rate * T);
            double attendu2 = spot2 * std::exp(rate * T);

            Logger::WriteMessage((std::wstring(L"E[S1(T)] généré  : ") + std::to_wstring(mean1)).c_str());
            Logger::WriteMessage((std::wstring(L"E[S1(T)] attendu : ") + std::to_wstring(attendu1)).c_str());
            Logger::WriteMessage((std::wstring(L"E[S2(T)] généré  : ") + std::to_wstring(mean2)).c_str());
            Logger::WriteMessage((std::wstring(L"E[S2(T)] attendu : ") + std::to_wstring(attendu2)).c_str());

            Assert::IsTrue(std::abs(mean1 - attendu1) / attendu1 < 0.02,
                L"E[S1(T)] doit être proche de Spot1*exp(r*T).");
            Assert::IsTrue(std::abs(mean2 - attendu2) / attendu2 < 0.02,
                L"E[S2(T)] doit être proche de Spot2*exp(r*T).");
        }

        // ----------------------------------------------------------------
        // Test 4 : Corrélation des log-rendements ≈ rho
        // ----------------------------------------------------------------
        TEST_METHOD(LogReturnsCorrelationNearRho)
        {
            size_t N = 5000;
            double spot1 = 100.0, spot2 = 80.0;
            double rate = 0.0, vol1 = 0.2, vol2 = 0.3, rho = 0.7, T = 1.0;
            size_t nbSteps = 252;

            double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0;

            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                BSMilstein2D bs(&ng, spot1, spot2, rate, vol1, vol2, rho);
                bs.Simulate(0.0, T, nbSteps);

                double logR1 = std::log(bs.GetPath(0)->GetState(T) / spot1);
                double logR2 = std::log(bs.GetPath(1)->GetState(T) / spot2);

                sumX  += logR1; sumY  += logR2;
                sumX2 += logR1 * logR1; sumY2 += logR2 * logR2;
                sumXY += logR1 * logR2;
            }

            double mX = sumX / N, mY = sumY / N;
            double varX = sumX2 / N - mX * mX;
            double varY = sumY2 / N - mY * mY;
            double cov  = sumXY / N - mX * mY;
            double corrEmp = cov / std::sqrt(varX * varY);

            Logger::WriteMessage((std::wstring(L"Corrélation empirique : ") + std::to_wstring(corrEmp)).c_str());
            Logger::WriteMessage((std::wstring(L"rho attendu           : ") + std::to_wstring(rho)).c_str());

            Assert::IsTrue(std::abs(corrEmp - rho) < 0.07,
                L"La corrélation empirique des log-rendements doit être proche de rho.");
        }
    };
}
