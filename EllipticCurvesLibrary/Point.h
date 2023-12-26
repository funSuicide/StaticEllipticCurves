#pragma once
#pragma warning( disable : 4127 )
#pragma warning( disable : 4146 )
#pragma warning( disable : 4244 )
#include <gmpxx.h>
#include <iostream> 

namespace EllipticCurvesLibrary
{
    class Point
    {
    private:
        mpz_class x, y;
        bool i;
    public:
        Point(mpz_class x, mpz_class y, bool i);
        mpz_class getX() const;
        mpz_class getY() const;
        bool onInfinity() const;
        Point operator=(Point p);
        bool operator==(const Point point) const;
       
        
        friend std::ostream& operator << (std::ostream& os, const Point& point);
      
    };

   

}