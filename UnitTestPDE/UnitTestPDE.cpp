#include "pch.h"
#include "CppUnitTest.h"
#include "../PDE/CallTerminalCondition.h"
#include "../PDE/PutTerminalCondition.h"
#include "../PDE/NullFunction.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestPDE
{
	TEST_CLASS(TerminalConditionTests)
	{
	public:
		
		TEST_METHOD(TestCallTerminalCondition)
		{
			CallTerminalCondition call(100.0);
			Assert::AreEqual(10.0, call(110.0), 1e-6);
			Assert::AreEqual(0.0, call(100.0), 1e-6);
			Assert::AreEqual(0.0, call(90.0), 1e-6);
		}

		TEST_METHOD(TestPutTerminalCondition)
		{
			PutTerminalCondition put(100.0);
			Assert::AreEqual(0.0, put(110.0), 1e-6);
			Assert::AreEqual(0.0, put(100.0), 1e-6);
			Assert::AreEqual(10.0, put(90.0), 1e-6);
		}

		TEST_METHOD(TestNullFunction)
		{
			NullFunction nullFunc;
			Assert::AreEqual(0.0, nullFunc(100.0, 1.0), 1e-6);
			Assert::AreEqual(0.0, nullFunc(-50.0, 0.5), 1e-6);
		}
	};
}
