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
	mpz_class x;
	mpz_class x2;
	mpz_class y;
	mpz_class y2;
	mpz_class n;
	std::string sForP = "499557256977219506380281345868505850806536870536066181139024969386512032771568813350165621784290905066210190413592751738596447503025705949107751116442344223236577235536504773278112648226610301283034800322975737256256519954201734040414098454787413957718966891611091832936397579180823123655431204786287";
	std::string sForA = "123123679126794129398123981293123821939129321391293120123123123123";
	std::string sForB = "1231231233284325982985987329749327432491283213812329839812932";
	std::string sX = "65536843059317579432679220396213332194376581011808234385727417169072394377738253318316784855818713245847263721284538469553508653876078761019991444556191420941978737598617327234255430558343815720203462123827413124075320030865091518906606073633044925237906801112803866248332473746948689125617050600609";
	std::string sX2 = "9796373221671279900825788503809012723532997697622734107257883110269968394493381895266796306412784898753425436538244326923250928";
	std::string sY = "4899072981123784531408949197300490762431742725924787564227633779352285858693210415539636076801067199854947992493207054116753763277608120043965346744875775487330626441906362612210553458453566405584366517604624896743865762779322612384290783376968180813062940363093347177276603914068861291024822026716";
	std::string sY2 = "25801540760717384198196093129334839815036182246776913720794951456209662946864473005489364001000966660950994855007127484625120425";
	std::string sN = "9839842989872364387643768876436980430498498435688432709843798";


	mpz_set_str(p.get_mpz_t(), sForP.c_str(), 10);
	mpz_set_str(a.get_mpz_t(), sForA.c_str(), 10);
	mpz_set_str(b.get_mpz_t(), sForB.c_str(), 10);

	mpz_set_str(x.get_mpz_t(), sX.c_str(), 10);
	mpz_set_str(x2.get_mpz_t(), sX2.c_str(), 10);
	mpz_set_str(y.get_mpz_t(), sY.c_str(), 10);
	mpz_set_str(y2.get_mpz_t(), sY2.c_str(), 10);
	mpz_set_str(n.get_mpz_t(), sN.c_str(), 10);


	std::cout << "Len P (bits): " << mpz_sizeinbase(p.get_mpz_t(), 2) << std::endl;
	EllipticCurvesLibrary::EllipticCurve E(a, b, p);
	std::cout << E << "\n" << std::endl;

	
	Point P(x, y, 0);
	Point Q = E.searchPoint();
	std::cout << "P: " << P << "\n" << std::endl;
	std::cout << "Q: " << Q << "\n" << std::endl;

	if (E.checkPoint(P)) {
		std::cout << "Point P is attached to curve" << "\n" << std::endl;
	}
	else
	{
		std::cout << "Point P is not attached to curve" << "\n" << std::endl;
	}

	if (E.checkPoint(Q)) {
		std::cout << "Point Q is attached to curve" << "\n" << std::endl;
	}
	else
	{
		std::cout << "Point Q is not attached to curve" << "\n" << std::endl;
	}

	Point R1 = E.sumTwoPoints(P, Q);

	std::cout << "R1: " << R1 << "\n" << std::endl;
	
	Point R2 = E.doubleAndAdd(P, n);
	Point R3 = E.ternaryMethod(P, n);

	std::cout << "R1: " << R1 << "\n" << std::endl;
	std::cout << "R2: " << R1 << "\n" << std::endl;
}
