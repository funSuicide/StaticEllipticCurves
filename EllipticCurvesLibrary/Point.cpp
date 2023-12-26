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
        return (x == point.x && i == point.i);
    }

    std::ostream& operator << (std::ostream& os, const Point& point)
    {
        std::string infinity = "normal";
        if (point.onInfinity())
        {
            infinity = "inf";
        }
        std::string outputString = "(" + point.x.get_str() + ", " + point.y.get_str() + ", " + infinity + ")";
        std::cout << outputString;
        return std::cout;
    }

}