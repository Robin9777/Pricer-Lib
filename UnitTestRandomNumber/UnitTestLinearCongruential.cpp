#include "pch.h"
#include "CppUnitTest.h"
#include <string>


#include "../RandomGenerator/LinearCongruential.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestRandomNumber
{

    TEST_CLASS(LinearCongruentialTests)
    {
    public:

        TEST_METHOD(FirstValueCorrect)
        {
            // Arrange : création du générateur avec des paramètres connus
            size_t seed = 1, a = 5, c = 3, m = 16;
            LinearCongruential gen(seed, a, c, m);

            // Act : on génère une valeur
            double valeur = gen.Generate();

            // Assert : on compare avec la valeur calculée à la main
            // X_1 = (5*1 + 3) % 16 = 8   =>   8 / 16 = 0.5
            double attendu = 8.0 / 16.0;
            Assert::AreEqual(attendu, valeur, 1e-12,
                L"La première valeur générée ne correspond pas à la formule X1 = (a*seed + c) % m.");

            Logger::WriteMessage((std::wstring(L"number generated : ") + std::to_wstring(valeur)).c_str());
            Logger::WriteMessage((std::wstring(L"expected number : ") + std::to_wstring(attendu)).c_str());
        }

        TEST_METHOD(SecondeValueCorrect)
        {
            size_t seed = 1, a = 5, c = 3, m = 16;
            LinearCongruential gen(seed, a, c, m);

            gen.Generate();               // X1 = 8  (on "consomme" la première)
            double valeur = gen.Generate(); // X2 = (5*8 + 3) % 16 = 11

            double attendu = 11.0 / 16.0;
            Assert::AreEqual(attendu, valeur, 1e-12,
                L"La deuxième valeur ne correspond pas à la récurrence chainée.");
        }

        TEST_METHOD(ValueIn01)
        {
            size_t seed = 42, a = 1664525, c = 1013904223, m = 4294967296ULL;
            LinearCongruential gen(seed, a, c, m);

            int N = 1000;
            for (int i = 0; i < N; ++i)
            {
                double v = gen.Generate();
                Assert::IsTrue(v >= 0.0 && v < 1.0,
                    L"Une valeur générée est en dehors de l'intervalle [0, 1).");
            }
        }


        TEST_METHOD(SameSeedSameSequence)
        {
            size_t seed = 7, a = 5, c = 3, m = 16;
            LinearCongruential gen1(seed, a, c, m);
            LinearCongruential gen2(seed, a, c, m);

            for (int i = 0; i < 10; ++i)
            {
                Assert::AreEqual(gen1.Generate(), gen2.Generate(), 1e-12,
                    L"Deux générateurs avec la même graine produisent des valeurs différentes.");
            }
        }

        TEST_METHOD(GetModulusReturnValue)
        {
            size_t m = 16;
            LinearCongruential gen(1, 5, 3, m);

            Assert::AreEqual(m, gen.GetModulus(),
                L"GetModulus() ne retourne pas la valeur passée au constructeur.");
        }
    };
}

