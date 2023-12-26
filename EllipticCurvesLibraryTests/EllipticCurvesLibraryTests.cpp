#include "CppUnitTest.h"
#include "../EllipticCurvesLibrary/EllipticCurve.h"
#include "../EllipticCurvesLibrary/EllipticCurve.cpp"
#include <string>

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
using namespace EllipticCurvesLibrary;

namespace Microsoft
{
	namespace VisualStudio
	{
		namespace CppUnitTestFramework
		{
			template<>
			static std::wstring ToString<Point>(const Point& point) {
				std::string coordX = point.getX().get_str();
				std::string coordY = point.getY().get_str();
				std::string output = "(" + coordX + ";" + coordY + ")";
				return std::wstring(output.begin(), output.end());
			}

			template<>
			static std::wstring ToString<std::vector<int>>(const std::vector<int>& v) {
				std::string output;
				for (auto i = v.begin(); i != v.end(); ++i) {
					switch (*i)
					{
					case -1:
						output += "-1, ";
						break;
					case 0:
						output += "0, ";
						break;
					case 1:
						output += "1, ";
						break;
					default:
						break;
					}
				}
				return std::wstring(output.begin(), output.end());
			}
		}
	}
}

namespace EllipticCurvesLibraryTests
{
	TEST_CLASS(EllipticCurvesLibraryTests)
	{
	public:
		TEST_METHOD(TestEEA)
		{
			mpz_class a(2);
			mpz_class p(107);
			std::vector<mpz_class> res = EEA(a, p);
			mpz_class actual = res[0];
			mpz_class expected(1);
			Assert::AreEqual(expected.get_str(10), actual.get_str(10));
		}
		TEST_METHOD(TestReversePoint)
		{
			EllipticCurve E(1, 3, 107);
			Point P(6, 92, 0);
			Point actual = E.reversePoint(P);
			Point expected(6, 15, 0);
			Assert::AreEqual(actual, expected);
		}

		TEST_METHOD(TestSumPoints)
		{
			EllipticCurve E(1, 3, 107);
			Point P(6, 92, 0);
			Point Q(2, 21, 0);
			Point actual = E.sumTwoPoints(P, Q);
			Point expected(73, 83, 0);
			Assert::AreEqual(expected, actual);
		}

		TEST_METHOD(TestDoubleAndAdd)
		{
			EllipticCurve E(1, 3, 107);
			Point P(6, 92, 0);
			mpz_class n(50);
			Point actual = E.doubleAndAdd(P, n);
			Point expected(106, 1, 0);
			Assert::AreEqual(expected, actual);
		}

		TEST_METHOD(testTernaryDecomposition)
		{
			mpz_class test = 57188;
			std::vector<int> expected = {1, 0, 0, -1, 0, 0, 0, 0, 0, -1, 0, -1, 0, 0, 1, 0, 0};
			std::vector<int> actual = ternaryDecomposition(test);
			Assert::AreEqual(expected, actual);
		}

		TEST_METHOD(testTernaryMehtod)
		{
			EllipticCurve E(-2, 7, 19);
			Point P(1, 5, 0);
			mpz_class n(15);
			Point actual = E.ternaryMethod(P, n);
			Point expected(106, 1, 0);
			Assert::AreEqual(expected, actual);
		}

		TEST_METHOD(testCurveCount)
		{

		}
	};
}
