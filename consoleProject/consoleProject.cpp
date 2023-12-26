#include "EllipticCurve.h"
#include <gmpxx.h>
#include <iostream>
#include <set>

using namespace EllipticCurvesLibrary;

int main()
{
	mpz_class p;
	mpz_class a;
	mpz_class b;
	mpz_class x1;
	mpz_class x2;
	mpz_class y1;
	mpz_class y2;
	mpz_class n;
	std::string sForP = "29933500047097369067937316004196040052873624302205671175238976482131293332251605769347243889197845099882939039929103427802794111";
	std::string sForA = "123123679126794129398123981293123821939129321391293120123123123123";
	std::string sForB = "1231231233284325982985987329749327432491283213812329839812932";
	std::string sX1 = "29244634081723867984438788561326585212900803192359039624106664878889385752078911721652892386416297514704264622266882233069590590";
	std::string sX2 = "9796373221671279900825788503809012723532997697622734107257883110269968394493381895266796306412784898753425436538244326923250928";
	std::string sY1 = "3103463007522824575400885612643449683050049481367289110460009970054903906870734032856822506722779226673508674520703832292023951";
	std::string sY2 = "25801540760717384198196093129334839815036182246776913720794951456209662946864473005489364001000966660950994855007127484625120425";
	std::string sN = "9839842989872364387643768876436980430498498435688432709843798";

	mpz_set_str(p.get_mpz_t(), sForP.c_str(), 10);
	mpz_set_str(a.get_mpz_t(), sForA.c_str(), 10);
	mpz_set_str(b.get_mpz_t(), sForB.c_str(), 10);

	mpz_set_str(x1.get_mpz_t(), sX1.c_str(), 10);
	mpz_set_str(x2.get_mpz_t(), sX2.c_str(), 10);
	mpz_set_str(y1.get_mpz_t(), sY1.c_str(), 10);
	mpz_set_str(y2.get_mpz_t(), sY2.c_str(), 10);
	mpz_set_str(n.get_mpz_t(), sN.c_str(), 10);
	
	EllipticCurvesLibrary::EllipticCurve E(a, b, p);
	Point P(x1, y1, 0);
	Point Q(x2, y2, 0);

	//Point R1 = E.doubleAndAdd(P, n);
	//Point R2 = E.ternaryMethod(P, n);

	//std::cout << R1.getX().get_str() << ", " << R1.getY().get_str() << std::endl;
	//std::cout << R2.getX().get_str() << ", " << R2.getY().get_str() << std::endl;
	mpz_class RES = E.curveCount();
	std::cout << "Order: " << RES.get_str() << std::endl;

	//std::set<mpz_class> S;
	//S.insert(b);
}
