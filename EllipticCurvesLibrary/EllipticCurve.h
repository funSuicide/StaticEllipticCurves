#pragma once
#include <math.h>
#include "Point.h"

namespace EllipticCurvesLibrary
{
    class EllipticCurve
    {
    private:
        mpz_class a, b, p;
    public:
        EllipticCurve() {};
        EllipticCurve(mpz_class a, mpz_class b, mpz_class p);
        Point sumTwoPoints(const Point& first, const Point& second) const;
        Point reversePoint(const Point& point) const;
        Point doubleAndAdd(const Point point, mpz_class n) const;
        Point ternaryMethod(const Point point, mpz_class n) const;
        Point searchPoint() const;
        bool checkPoint(Point& point) const;
        mpz_class curveOrder() const;
        mpz_class getA() const;
        mpz_class getB() const;

        friend std::ostream& operator << (std::ostream& os, const EllipticCurve& );
    };
}