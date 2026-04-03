#include "pch.h"
#include "CppUnitTest.h"
#include <string>
#include <cmath>
#include <vector>

#include "../SDE/BrownianD1.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/Normal.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestSDE
{
    TEST_CLASS(Brownian1DTests)
    {
    public:

        // ----------------------------------------------------------------
        // Helper : générateur N(0,1) via LinearCongruential + Box-Muller
        // ----------------------------------------------------------------
        static Normal MakeNormalGen(size_t seed = 42)
        {
            // Paramètres LCG classiques (Numerical Recipes)
            static LinearCongruential lc(seed, 1664525, 1013904223, 4294967296ULL);
            static Normal ng(0.0, 1.0, lc);
            return ng;  // retourne par valeur car il n'y a pas d'état interne partagé critique ici
        }

        // ----------------------------------------------------------------
        // Test 1 : W(0) = 0  (valeur initiale)
        // ----------------------------------------------------------------
        TEST_METHOD(InitialValueIsZero)
        {
            LinearCongruential lc(42, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            Brownian1D W(&ng);
            W.Simulate(0.0, 1.0, 100);

            double w0 = W.GetPath(0)->GetState(0.0);

            Logger::WriteMessage((std::wstring(L"W(0) génered : ") + std::to_wstring(w0)).c_str());
            Logger::WriteMessage((std::wstring(L"W(0) expected : ") + std::to_wstring(0.0)).c_str());

            Assert::AreEqual(0.0, w0, 1e-12, L"W(0) doit être 0.");
        }

        // ----------------------------------------------------------------
        // Test 2 : Le nombre de valeurs dans le path = nbSteps + 1
        // ----------------------------------------------------------------
        TEST_METHOD(PathSizeIsNbStepsPlusOne)
        {
            LinearCongruential lc(1, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            size_t nbSteps = 50;
            Brownian1D W(&ng);
            W.Simulate(0.0, 1.0, nbSteps);

            std::vector<double> vals = W.GetPath(0)->GetAllValues();
            size_t valeur = vals.size();
            size_t attendu = nbSteps + 1;

            Logger::WriteMessage((std::wstring(L"Taille du path générée : ") + std::to_wstring(valeur)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille du path attendue : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, valeur, L"Le path doit contenir nbSteps + 1 valeurs.");
        }

        // ----------------------------------------------------------------
        // Test 3 : Propriété de martingale : E[W(T)] ≈ 0
        //          On simule N chemins et on moyenne W(T).
        // ----------------------------------------------------------------
        TEST_METHOD(MartingaleProperty_MeanIsZero)
        {
            size_t N = 5000;
            double T = 1.0;
            size_t nbSteps = 100;

            double sum = 0.0;
            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                Brownian1D W(&ng);
                W.Simulate(0.0, T, nbSteps);
                sum += W.GetPath(0)->GetState(T);
            }

            double meanWT = sum / N;

            Logger::WriteMessage((std::wstring(L"E[W(T)] généré  : ") + std::to_wstring(meanWT)).c_str());
            Logger::WriteMessage((std::wstring(L"E[W(T)] attendu : ") + std::to_wstring(0.0)).c_str());

            // Tolérance : 3 * sigma / sqrt(N)  avec sigma=1 => ~0.042
            Assert::IsTrue(std::abs(meanWT) < 0.05,
                L"E[W(T)] devrait être proche de 0 (propriété de martingale).");
        }

        // ----------------------------------------------------------------
        // Test 4 : Variance de W(T) ≈ T  (propriété fondamentale du BM)
        // ----------------------------------------------------------------
        TEST_METHOD(VarianceEqualsT)
        {
            size_t N = 5000;
            double T = 2.0;
            size_t nbSteps = 200;

            double sum = 0.0, sumSq = 0.0;
            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                Brownian1D W(&ng);
                W.Simulate(0.0, T, nbSteps);
                double wT = W.GetPath(0)->GetState(T);
                sum += wT;
                sumSq += wT * wT;
            }

            double mean = sum / N;
            double varWT = sumSq / N - mean * mean;

            Logger::WriteMessage((std::wstring(L"Var[W(T)] générée  : ") + std::to_wstring(varWT)).c_str());
            Logger::WriteMessage((std::wstring(L"Var[W(T)] attendue : ") + std::to_wstring(T)).c_str());

            Assert::IsTrue(std::abs(varWT - T) < 0.15,
                L"Var[W(T)] devrait être proche de T.");
        }

        // ----------------------------------------------------------------
        // Test 5 : ré-simulation nette (Simulate efface le path précédent)
        // ----------------------------------------------------------------
        TEST_METHOD(ResimulationClearsPreviousPath)
        {
            LinearCongruential lc(7, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            Brownian1D W(&ng);
            W.Simulate(0.0, 1.0, 50);
            W.Simulate(0.0, 1.0, 30);   // ré-simulation avec nbSteps plus petit

            std::vector<double> vals = W.GetPath(0)->GetAllValues();
            size_t valeur = vals.size();
            size_t attendu = 30 + 1;

            Logger::WriteMessage((std::wstring(L"Taille après ré-simulation : ") + std::to_wstring(valeur)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille attendue            : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, valeur,
                L"Après ré-simulation le path doit avoir la nouvelle taille.");
        }
    };
}
