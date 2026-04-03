#include "pch.h"
#include "CppUnitTest.h"
#include <string>
#include <cmath>
#include <vector>

#include "../SDE/BrownianND.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/Normal.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestSDE
{
    TEST_CLASS(BrownianNDTests)
    {
    public:

        // ----------------------------------------------------------------
        // Test 1 : Chaque dimension démarre à 0
        // ----------------------------------------------------------------
        TEST_METHOD(InitialValuesAreZero)
        {
            LinearCongruential lc(42, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            std::vector<std::vector<double>> corrId = { {1.0, 0.0}, {0.0, 1.0} };
            BrownianND W(&ng, 2, &corrId);
            W.Simulate(0.0, 1.0, 100);

            double w0_dim0 = W.GetPath(0)->GetState(0.0);
            double w0_dim1 = W.GetPath(1)->GetState(0.0);

            Logger::WriteMessage((std::wstring(L"W[0](0) généré  : ") + std::to_wstring(w0_dim0)).c_str());
            Logger::WriteMessage((std::wstring(L"W[1](0) généré  : ") + std::to_wstring(w0_dim1)).c_str());

            Assert::AreEqual(0.0, w0_dim0, 1e-12, L"W[0](0) doit être 0.");
            Assert::AreEqual(0.0, w0_dim1, 1e-12, L"W[1](0) doit être 0.");
        }

        // ----------------------------------------------------------------
        // Test 2 : Nombre de valeurs = nbSteps + 1 pour chaque dimension
        // ----------------------------------------------------------------
        TEST_METHOD(PathSizeCorrectForEachDimension)
        {
            LinearCongruential lc(1, 1664525, 1013904223, 4294967296ULL);
            Normal ng(0.0, 1.0, lc);

            std::vector<std::vector<double>> corrId = { {1.0, 0.0}, {0.0, 1.0} };
            size_t nbSteps = 40;
            BrownianND W(&ng, 2, &corrId);
            W.Simulate(0.0, 1.0, nbSteps);

            size_t n0 = W.GetPath(0)->GetAllValues().size();
            size_t n1 = W.GetPath(1)->GetAllValues().size();
            size_t attendu = nbSteps + 1;

            Logger::WriteMessage((std::wstring(L"Taille dim 0 : ") + std::to_wstring(n0)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille dim 1 : ") + std::to_wstring(n1)).c_str());
            Logger::WriteMessage((std::wstring(L"Taille attendue : ") + std::to_wstring(attendu)).c_str());

            Assert::AreEqual(attendu, n0, L"Dim 0 : path size doit être nbSteps+1.");
            Assert::AreEqual(attendu, n1, L"Dim 1 : path size doit être nbSteps+1.");
        }

        // ----------------------------------------------------------------
        // Test 3 : Matrice identité => processus indépendants (corrélation ≈ 0)
        // ----------------------------------------------------------------
        TEST_METHOD(UncorrelatedBrownians_CorrelationNearZero)
        {
            size_t N = 5000;
            double T = 1.0;
            size_t nbSteps = 100;

            double sumXY = 0.0, sumX = 0.0, sumY = 0.0;
            double sumX2 = 0.0, sumY2 = 0.0;

            std::vector<std::vector<double>> corrId = { {1.0, 0.0}, {0.0, 1.0} };

            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                BrownianND W(&ng, 2, &corrId);
                W.Simulate(0.0, T, nbSteps);

                double x = W.GetPath(0)->GetState(T);
                double y = W.GetPath(1)->GetState(T);
                sumX  += x;  sumY  += y;
                sumX2 += x*x; sumY2 += y*y;
                sumXY += x * y;
            }

            double meanX = sumX / N, meanY = sumY / N;
            double varX  = sumX2 / N - meanX * meanX;
            double varY  = sumY2 / N - meanY * meanY;
            double cov   = sumXY / N - meanX * meanY;
            double corr  = cov / std::sqrt(varX * varY);

            Logger::WriteMessage((std::wstring(L"Corrélation empirique : ") + std::to_wstring(corr)).c_str());
            Logger::WriteMessage((std::wstring(L"Corrélation attendue  : ") + std::to_wstring(0.0)).c_str());

            Assert::IsTrue(std::abs(corr) < 0.08,
                L"Les browniens indépendants doivent avoir une corrélation proche de 0.");
        }

        // ----------------------------------------------------------------
        // Test 4 : Matrice avec rho=0.8 => corrélation empirique ≈ 0.8
        // ----------------------------------------------------------------
        TEST_METHOD(CorrelatedBrownians_CorrelationNearRho)
        {
            size_t N = 5000;
            double T = 1.0;
            size_t nbSteps = 100;
            double rho = 0.8;

            std::vector<std::vector<double>> corrMat = { {1.0, rho}, {rho, 1.0} };

            double sumXY = 0.0, sumX = 0.0, sumY = 0.0;
            double sumX2 = 0.0, sumY2 = 0.0;

            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                BrownianND W(&ng, 2, &corrMat);
                W.Simulate(0.0, T, nbSteps);

                double x = W.GetPath(0)->GetState(T);
                double y = W.GetPath(1)->GetState(T);
                sumX  += x;  sumY  += y;
                sumX2 += x*x; sumY2 += y*y;
                sumXY += x * y;
            }

            double meanX = sumX / N, meanY = sumY / N;
            double varX  = sumX2 / N - meanX * meanX;
            double varY  = sumY2 / N - meanY * meanY;
            double cov   = sumXY / N - meanX * meanY;
            double corrEmpirique = cov / std::sqrt(varX * varY);

            Logger::WriteMessage((std::wstring(L"Corrélation empirique : ") + std::to_wstring(corrEmpirique)).c_str());
            Logger::WriteMessage((std::wstring(L"Corrélation attendue  : ") + std::to_wstring(rho)).c_str());

            Assert::IsTrue(std::abs(corrEmpirique - rho) < 0.05,
                L"La corrélation empirique doit être proche de rho=0.8.");
        }

        // ----------------------------------------------------------------
        // Test 5 : Variance de chaque dimension ≈ T (processus brownien standard)
        // ----------------------------------------------------------------
        TEST_METHOD(EachDimensionVarianceEqualsT)
        {
            size_t N = 5000;
            double T = 1.5;
            size_t nbSteps = 150;

            std::vector<std::vector<double>> corrId = { {1.0, 0.0}, {0.0, 1.0} };

            double sum0 = 0.0, sumSq0 = 0.0;
            double sum1 = 0.0, sumSq1 = 0.0;

            for (size_t s = 0; s < N; ++s)
            {
                LinearCongruential lc(s + 1, 1664525, 1013904223, 4294967296ULL);
                Normal ng(0.0, 1.0, lc);
                BrownianND W(&ng, 2, &corrId);
                W.Simulate(0.0, T, nbSteps);

                double w0 = W.GetPath(0)->GetState(T);
                double w1 = W.GetPath(1)->GetState(T);
                sum0   += w0; sumSq0 += w0 * w0;
                sum1   += w1; sumSq1 += w1 * w1;
            }

            double var0 = sumSq0 / N - (sum0 / N) * (sum0 / N);
            double var1 = sumSq1 / N - (sum1 / N) * (sum1 / N);

            Logger::WriteMessage((std::wstring(L"Var dim 0 : ") + std::to_wstring(var0)).c_str());
            Logger::WriteMessage((std::wstring(L"Var dim 1 : ") + std::to_wstring(var1)).c_str());
            Logger::WriteMessage((std::wstring(L"T attendu : ") + std::to_wstring(T)).c_str());

            Assert::IsTrue(std::abs(var0 - T) < 0.15, L"Var[W0(T)] doit être proche de T.");
            Assert::IsTrue(std::abs(var1 - T) < 0.15, L"Var[W1(T)] doit être proche de T.");
        }
    };
}
