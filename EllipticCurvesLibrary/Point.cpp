#include "Point.h"

namespace EllipticCurvesLibrary
{
    Point::Point(mpz_class x, mpz_class y, bool i)
    {
        this->x = x;
        this->y = y;
        this->i = i;
    }

    mpz_class Point::getX() const
    {
        return this->x;
    }

    mpz_class Point::getY() const
    {
        return this->y;
    }

    bool Point::onInfinity() const
    {
        return i;
    }

    Point Point::operator=(Point p)
    {
        x = p.getX();
        y = p.getY();
        i = p.onInfinity();
        return *this; 
    }

    bool Point::operator==(const Point point) const
    {
        return (x == point.getX() && i == point.onInfinity());
    }
}