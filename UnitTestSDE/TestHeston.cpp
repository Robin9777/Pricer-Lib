#include "pch.h"
#include "CppUnitTest.h"
#include <string>
#include <cmath>
#include <vector>

#include "../SDE/Heston.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/Normal.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestSDE
{
    TEST_CLASS(HestonTests)
    {
    public:
        // Paramètres Heston typiques utilisés dans plusieurs tests ci-dessous :
        //   spot=100, v0=0.04, mu=0.05, theta=0.04, kappa=2.0, sigma=0.3, rho=-0.7
        //   => v_bar = theta = 0.04, Feller: 2*kappa*theta = 0.16 > sigma^2 = 0.09 (OK)

        // ----------------------------------------------------------------
        // Test 1 : S(0) = spot, V(0) = v0  (valeurs initiales dans les paths)
        // ----------------------------------------------------------------
        TEST_METHOD(InitialValuesAreCorrect)
        {
            LinearCongruential lc(42, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            double spot = 100.0, v0 = 0.04;
            Heston heston(&ng, spot, v0, 0.05, 0.04, 2.0, 0.3, -0.7);
            heston.Simulate(0.0, 1.0, 200);

            double s0 = heston.GetPath(0)->GetState(0.0);
            double v0_path = heston.GetPath(1)->GetState(0.0);  // path[1] = variance

            Logger::WriteMessage((std::wstring(L"S(0) généré  : ") + std::to_wstring(s0)).c_str());
            Logger::WriteMessage((std::wstring(L"S(0) attendu : ") + std::to_wstring(spot)).c_str());
            Logger::WriteMessage((std::wstring(L"V(0) généré  : ") + std::to_wstring(v0_path)).c_str());
            Logger::WriteMessage((std::wstring(L"V(0) attendu : ") + std::to_wstring(v0)).c_str());

            Assert::AreEqual(spot, s0,     1e-10, L"S(0) doit être le spot initial.");
            Assert::AreEqual(v0,  v0_path, 1e-10, L"V(0) doit être la variance initiale.");
        }

        // ----------------------------------------------------------------
        // Test 2 : Les deux paths ont la taille nbSteps + 1
        // ----------------------------------------------------------------
        TEST_METHOD(BothPathSizesAreCorrect)
        {
            LinearCongruential lc(1, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            size_t nbSteps = 120;
            Heston heston(&ng, 100.0, 0.04, 0.05, 0.04, 2.0, 0.3, -0.7);
            heston.Simulate(0.0, 1.0, nbSteps);

            size_t nS = heston.GetPath(0)->GetAllValues().size();
            size_t nV = heston.GetPath(1)->GetAllValues().size();
            size_t attendu = nbSteps + 1;

            Logger::WriteMessage((std::wstring(L"Taille path S : ") + std::to_wstring(nS)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille path V : ") + std::to_wstring(nV)).c_str());
            Logger::WriteMessage((std::wstring(L"Attendu       : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, nS, L"Path S : taille doit être nbSteps+1.");
            Assert::AreEqual(attendu, nV, L"Path V : taille doit être nbSteps+1.");
        }

        // ----------------------------------------------------------------
        // Test 3 : E[S(T)] ≈ S0 * exp(mu*T)  (drift risk-neutral)
        // ----------------------------------------------------------------
        TEST_METHOD(MeanPriceConvergesToForwardPrice)
        {
            size_t N = 5000;
            double spot = 100.0, mu = 0.05, T = 1.0;
            size_t nbSteps = 250;

            double sum = 0.0;
            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                Heston heston(&ng, spot, 0.04, mu, 0.04, 2.0, 0.3, -0.7);
                heston.Simulate(0.0, T, nbSteps);
                sum += heston.GetPath(0)->GetState(T);
            }

            double mean    = sum / N;
            double attendu = spot * std::exp(mu * T);

            Logger::WriteMessage((std::wstring(L"E[S(T)] généré  : ") + std::to_wstring(mean)).c_str());
            Logger::WriteMessage((std::wstring(L"E[S(T)] attendu : ") + std::to_wstring(attendu)).c_str());

            Assert::IsTrue(std::abs(mean - attendu) / attendu < 0.03,
                L"E[S(T)] doit être proche de S0*exp(mu*T).");
        }

        // ----------------------------------------------------------------
        // Test 4 : Mean-reversion de la variance : E[V(T)] ≈ theta + (v0-theta)*exp(-kappa*T)
        // ----------------------------------------------------------------
        TEST_METHOD(VarianceMeanReversion)
        {
            size_t N = 5000;
            double v0 = 0.09, theta = 0.04, kappa = 3.0, sigma = 0.3;
            double T = 1.0;
            size_t nbSteps = 250;

            double sum = 0.0;
            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                // v0 > theta => on part "trop haut" et la variance doit revenir vers theta
                Heston heston(&ng, 100.0, v0, 0.0, theta, kappa, sigma, -0.5);
                heston.Simulate(0.0, T, nbSteps);
                sum += heston.GetPath(1)->GetState(T);
            }

            double meanV   = sum / N;
            double attendu = theta + (v0 - theta) * std::exp(-kappa * T);

            Logger::WriteMessage((std::wstring(L"E[V(T)] généré  : ") + std::to_wstring(meanV)).c_str());
            Logger::WriteMessage((std::wstring(L"E[V(T)] attendu : ") + std::to_wstring(attendu)).c_str());

            // Tolérance plus large car le CIR a une variance non nulle
            Assert::IsTrue(std::abs(meanV - attendu) < 0.01,
                L"E[V(T)] doit converger vers la formule analytique CIR.");
        }

        // ----------------------------------------------------------------
        // Test 5 : S(t) > 0  (le prix ne devient pas négatif)
        // ----------------------------------------------------------------
        TEST_METHOD(PriceRemainsPositive)
        {
            size_t N = 1000;
            double T = 1.0;
            size_t nbSteps = 252;

            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                Heston heston(&ng, 100.0, 0.04, 0.05, 0.04, 2.0, 0.3, -0.7);
                heston.Simulate(0.0, T, nbSteps);

                std::vector<double> vals = heston.GetPath(0)->GetAllValues();
                for (double v : vals)
                {
                    Assert::IsTrue(v > 0.0, L"S(t) ne doit jamais être négatif.");
                }
            }
        }

        // ----------------------------------------------------------------
        // Test 6 : Ré-simulation — les anciens paths sont effacés
        // ----------------------------------------------------------------
        TEST_METHOD(ResimulationClearsPreviousPath)
        {
            LinearCongruential lc(5, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            Heston heston(&ng, 100.0, 0.04, 0.05, 0.04, 2.0, 0.3, -0.7);
            heston.Simulate(0.0, 1.0, 200);
            heston.Simulate(0.0, 1.0, 50);   // ré-simulation avec nbSteps plus petit

            size_t valeur  = heston.GetPath(0)->GetAllValues().size();
            size_t attendu = 50 + 1;

            Logger::WriteMessage((std::wstring(L"Taille après ré-simulation : ") + std::to_wstring(valeur)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille attendue            : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, valeur,
                L"Après ré-simulation, le path doit avoir la nouvelle taille.");
        }

        // ----------------------------------------------------------------
        // Test 7 : sigma=0 => variance constante (CIR sans diffusion)
        //          V(t) = v0 + kappa*(theta-v0)*t  (Euler déterministe)
        //          On teste que V reste proche de theta après convergence.
        // ----------------------------------------------------------------
        TEST_METHOD(ZeroSigma_VarianceConvergesToTheta)
        {
            LinearCongruential lc(99, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            double v0 = 0.04, theta = 0.04, kappa = 2.0, sigma = 0.0, T = 1.0;
            size_t nbSteps = 1000;

            Heston heston(&ng, 100.0, v0, 0.0, theta, kappa, sigma, 0.0);
            heston.Simulate(0.0, T, nbSteps);

            double vT     = heston.GetPath(1)->GetState(T);
            double attendu = theta;   // CIR sans bruit => V reste à theta si v0=theta

            Logger::WriteMessage((std::wstring(L"V(T) généré  : ") + std::to_wstring(vT)).c_str());
            Logger::WriteMessage((std::wstring(L"V(T) attendu : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, vT, 1e-8,
                L"Avec sigma=0 et v0=theta, V(T) doit rester exactement à theta.");
        }
    };
}
