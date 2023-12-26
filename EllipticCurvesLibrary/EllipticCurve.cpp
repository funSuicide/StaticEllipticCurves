#pragma once
#include "EllipticCurve.h"
#include <string>
#include <vector>
#include <sstream>
#include <set>
#include <iostream>

mpz_class moduleTransform(mpz_class value, mpz_class p)
{
    return (value % p + p) % p;
}

std::vector<mpz_class> EEA(mpz_class a, mpz_class b)
{
    std::vector<mpz_class> results(2);
    mpz_class s = 0;
    mpz_class oldS = 1;
    mpz_class r = b;
    mpz_class oldR = a;

    while (r != 0)
    {
        mpz_class quotient = oldR / r;
        mpz_class tmp = oldR;
        oldR = r;
        r = tmp - quotient * r;
        tmp = oldS;
        oldS = s;
        s = tmp - quotient * s;
    }

    results[0] = oldR;
    results[1] = oldS;

    return results;
}

mpz_class reverseValue(mpz_class value, mpz_class p)
{
    std::vector<mpz_class> resultsEEA = EEA(value, p);

    mpz_class gcd = resultsEEA[0];
    mpz_class x = resultsEEA[1];
    if (gcd != 1)
    {
        // exept
    }
    return moduleTransform(x, p);
}

std::pair<mpz_class, mpz_class> getTAndS(mpz_class value)
{
    mpz_class copyValue(value);
    mpz_class s(0);
    mpz_class t(0);

    while (copyValue % 2 == 0)
    {
        copyValue /= 2;
        s++;
    }
    t = copyValue;
    return std::pair<mpz_class, mpz_class>(t, s);
}

mpz_class squareModule(mpz_class p, mpz_class a)
{
    gmp_randstate_t state;
    mpz_class seed;
    seed = (int)time(NULL);
    gmp_randinit_mt(state);
    gmp_randseed(state, seed.get_mpz_t());

    mpz_class d;

    while (true)
    {
        mpz_urandomm(d.get_mpz_t(), state, mpz_class(p - 1).get_mpz_t());
        if (d > 1) {
            if (mpz_legendre(d.get_mpz_t(), p.get_mpz_t()) == -1) break;
        }
    }

    mpz_class A(0);
    mpz_class D(0);
    std::pair<mpz_class, mpz_class> tAndS = getTAndS(p - 1);
    __gmpz_powm(A.get_mpz_t(), a.get_mpz_t(), tAndS.first.get_mpz_t(), p.get_mpz_t());
    __gmpz_powm(D.get_mpz_t(), d.get_mpz_t(), tAndS.first.get_mpz_t(), p.get_mpz_t());
    mpz_class m = 0;
    for (int i = 1; i < tAndS.second; ++i)
    {
        mpz_class tmp;
        __gmpz_powm(tmp.get_mpz_t(), D.get_mpz_t(), m.get_mpz_t(), p.get_mpz_t());
        tmp = A * tmp;//moduleTransform(A * D, p);

        mpz_class tmp2;
        mpz_class tmp3 = tAndS.second - 1 - i; //first -> second
        __gmpz_powm(tmp2.get_mpz_t(), mpz_class(2).get_mpz_t(), tmp3.get_mpz_t(), p.get_mpz_t());
        __gmpz_powm(tmp.get_mpz_t(), tmp.get_mpz_t(), tmp2.get_mpz_t(), p.get_mpz_t());
        if (tmp == moduleTransform(-1, p))
        {
            mpz_class tmpM;
            __gmpz_powm(tmpM.get_mpz_t(), mpz_class(2).get_mpz_t(), mpz_class(i).get_mpz_t(), p.get_mpz_t());
            m += tmpM;
        }
    }

    mpz_class tmp4;
    mpz_class tmp4S = moduleTransform((tAndS.first + 1) * reverseValue(2, p), p);
    __gmpz_powm(tmp4.get_mpz_t(), a.get_mpz_t(), tmp4S.get_mpz_t(), p.get_mpz_t());
    mpz_class tmp5;
    mpz_class tmp5S = moduleTransform(m * reverseValue(2, p), p); // m / 2;
    __gmpz_powm(tmp5.get_mpz_t(), D.get_mpz_t(), tmp5S.get_mpz_t(), p.get_mpz_t());
    mpz_class x = moduleTransform(tmp4 * tmp5, p);
    return x;
}

std::vector<int> ternaryDecomposition(mpz_class& n)
{
    int k = mpz_sizeinbase(n.get_mpz_t(), 2);
    std::vector<int> resultDecomposition(k+1);
    mpz_class bit;
    int count = 0;
    bool flag = 0;

    for (int i = 0; i < k; ++i)
    {
        bit = (n >> i) & 1;
        resultDecomposition[k - i] = bit.get_d();
    }
    resultDecomposition[0] = 0;

    for (int i = 0; i < k; ++i)
    {
        if ((resultDecomposition[k - i] == 1) && (resultDecomposition[k - i - 1] == 1) && !flag)
        {
            resultDecomposition[k - i] = -1;
            flag = 1;
        }
        else if ((resultDecomposition[k - i] == 1) && (resultDecomposition[k - i - 1] == 1) && flag)
        {
            resultDecomposition[k - i] = 0;
        }
        else if ((resultDecomposition[k - i] == 1) && (resultDecomposition[k - i - 1] == 0) && flag)
        {
            resultDecomposition[k - i] = 0;
            resultDecomposition[k - i - 1] = 1;
            flag = false;
        }
    }
    return resultDecomposition;
}

mpz_class getNativeSum(const mpz_class& a, const mpz_class&b, const mpz_class& p)
{
    mpz_class sum = 0;
    for (mpz_class x = 0; x < p; ++x)
    {
        mpz_class tmp = x * x * x + x * x * a + b;
        sum += mpz_legendre(tmp.get_mpz_t(), p.get_mpz_t());
    }
    return sum;
}




namespace EllipticCurvesLibrary
{
    EllipticCurve::EllipticCurve(mpz_class a, mpz_class b, mpz_class p)
    {
        this->a = a;
        this->b = b;
        this->p = p;
    }

    Point EllipticCurve::sumTwoPoints(const Point& first, const Point& second) const
    {
        mpz_class x1 = moduleTransform(first.getX(), p);
        mpz_class x2 = moduleTransform(second.getX(), p);
        mpz_class y1 = moduleTransform(first.getY(), p);
        mpz_class y2 = moduleTransform(second.getY(), p);
        mpz_class lambda = 0;

        if ((x1 == x2) && (y1 == moduleTransform((-y2), p)))
        {
            return Point(0, 0, 1);
        }
        if (first.onInfinity())
        {
            return second;
           
        } 
        if (second.onInfinity())
        {
            return first;
        }

        if (first == second)
        {
            lambda = moduleTransform((3 * x1 * x1 + this->a) * reverseValue(moduleTransform(2*y1, p), p), p);
        }
        else
        {
            lambda = moduleTransform((y2 - y1) * reverseValue(moduleTransform(x2 - x1, p), p), p);
        }
       
        mpz_class x3 = moduleTransform(lambda * lambda - x1 - x2, p);
        mpz_class y3 = moduleTransform(lambda * (x1- x3) - y1, p);
        return Point(x3, y3, 0);
    }

    Point EllipticCurve::reversePoint(const Point& point) const
    {
        if (point.onInfinity())
        {
            return Point(0, 0, 1);
        }
        else
        {
            mpz_class x = moduleTransform(point.getX(), p);
            mpz_class y = moduleTransform(-point.getY(), p);
            return Point(x, y, 0);
        }
    }

    Point EllipticCurve::doubleAndAdd(const Point point, mpz_class n) const {
        Point T = point;
        mpz_class bit;
        int sizeN = mpz_sizeinbase(n.get_mpz_t(), 2);
        Point Q(0, 0, 1);
        for (int i = 0; i < sizeN; ++i)
        {
            bit = (n >> i) & 1;
            if (bit == 1)
            {
                Q = sumTwoPoints(Q, T);
            }
            T = sumTwoPoints(T, T);
        }
        return Q;
    }

    Point EllipticCurve::ternaryMethod(const Point point, mpz_class n) const {
        Point T = point;
        Point Q(0, 0, 1);
        std::vector<int> ternaryDecompose = ternaryDecomposition(n);
        for (int i = 0; i < ternaryDecompose.size(); ++i) {
            if (ternaryDecompose[ternaryDecompose.size() - i - 1] == 1) {
                Q = sumTwoPoints(Q, T);
            }
            else if (ternaryDecompose[ternaryDecompose.size() - i - 1] == -1)
            {
                Q = sumTwoPoints(Q, reversePoint(T));
            }
            T = sumTwoPoints(T, T);
        }
        return Q;
    }

    Point EllipticCurve::searchPoint() const
    {
        gmp_randstate_t state;
        mpz_class t;
        mpz_class seed;
        seed = (int)time(NULL);
        gmp_randinit_mt(state);
        gmp_randseed(state, seed.get_mpz_t());
        mpz_class x;

      
        while (true)
        {
            mpz_urandomm(x.get_mpz_t(), state, mpz_class(p-1).get_mpz_t());
            mpz_class tmp = (x * (x * x + a) + b);
            t = moduleTransform(tmp, this->p);
            if (mpz_legendre(t.get_mpz_t(), p.get_mpz_t()) != -1)
            {
                mpz_class first = x;
                mpz_class second = squareModule(p, t);
                return Point(first, second, 0);
            }
        }
    }

    bool EllipticCurve::checkPoint(Point& point) const
    {
        mpz_class x = point.getX();
        mpz_class y = point.getY();
        return (moduleTransform(y * y, p) == moduleTransform(x * x * x + a * x + b, p));
    }

    mpz_class EllipticCurve::getA() const
    {
        return this->a;
    }

    mpz_class EllipticCurve::getB() const
    {
        return this->b;
    }

    mpz_class EllipticCurve::curveCount() const
    {
        gmp_randstate_t state;
        mpz_class g;
        mpz_class seed;
        seed = (int)time(NULL);
        gmp_randinit_mt(state);
        gmp_randseed(state, seed.get_mpz_t());

        mpz_class W; // параметр для гигантских шагов
        mpz_class c;
        mpz_class d; //параметры искажения

        if (p <= 229)
        {
            return p + 1 + getNativeSum(a, b, p);
        }
        else
        {
            mpz_urandomm(g.get_mpz_t(), state, p.get_mpz_t());
            g = 823; ///
            while (mpz_legendre(g.get_mpz_t(), p.get_mpz_t()) != -1)
            {
                mpz_urandomm(g.get_mpz_t(), state, p.get_mpz_t());
            }
            mpf_class floatP;
            mpf_set_z(floatP.get_mpf_t(), p.get_mpz_t());
            __gmpf_set_default_prec(1000);
            floatP = std::pow(floatP.get_d(), 0.25);
            floatP *= sqrt(2);
            __gmpf_ceil(floatP.get_mpf_t(), floatP.get_mpf_t());
            W = floatP;

            std::cout << W.get_str() << std::endl;

            __gmpz_powm(c.get_mpz_t(), g.get_mpz_t(), mpz_class(2).get_mpz_t(), p.get_mpz_t());
            c = moduleTransform(c * a, p);

            __gmpz_powm(d.get_mpz_t(), g.get_mpz_t(), mpz_class(3).get_mpz_t(), p.get_mpz_t());
            d = moduleTransform(d * b, p);

            mpz_class x;
            int sigma;
            while (true)
            {
                mpz_urandomm(x.get_mpz_t(), state, mpz_class(p-1).get_mpz_t());

                x = 400; /// 
                mpz_class tmp = x * x * x + x * a + b;
                sigma = mpz_legendre(tmp.get_mpz_t(), p.get_mpz_t());

                if (sigma != 0)
                {
                    EllipticCurve E = *this;
                    if (sigma  != 1)
                    {
                        E.a = c;
                        E.b = d;
                        x = moduleTransform(g * x, p);
                    }
                    mpz_class tmp1 = x * x * x + x * E.getA() + E.getB();
                    tmp1 = moduleTransform(tmp1, p);
                    mpz_class yForP = squareModule(p, tmp1);
                    Point P(x, 307, 0);

                    std::vector<mpz_class> S;
                    std::set<mpz_class> A;
                    std::set<mpz_class> B;
                    std::vector<mpz_class> tmpA;
                    std::vector<mpz_class> tmpB;

                    //std::cout << "List A: " << std::endl;
                    for (mpz_class beta = 0; beta <= W-1; ++beta)
                    {
                        mpz_class tmp = (p + 1 + beta);
                        Point setElement = E.doubleAndAdd(P, tmp);
                        mpz_class tmpElement = setElement.getX();
                       // std::cout << tmpElement.get_str() << std::endl;
                        A.insert(tmpElement);
                        tmpA.push_back(tmpElement);
                    }

                   // std::cout << "List B: " << std::endl;
                    for (mpz_class y = 0; y <= W; ++y)
                    {
                        mpz_class tmp = y * W;
                        Point setElement = E.doubleAndAdd(P, tmp);
                        mpz_class tmpElement = setElement.getX();
                       // std::cout << tmpElement.get_str() << std::endl;
                        B.insert(tmpElement);
                        tmpB.push_back(tmpElement);
                    }

                  
                    std::set_intersection(A.begin(), A.end(), B.begin(), B.end(), std::back_inserter(S));
                    
                    if (S.size() == 1)
                    {
                        Point O(0, 0, 1);
                        mpz_class s = S[0];

                        int indexA = 0;
                        int indexB = 0;

                        for (size_t i = 0; i < tmpA.size(); ++i)
                        {
                            if (tmpA[i] == s) indexA = i;
                        }
                        for (size_t i = 0; i < tmpB.size(); ++i)
                        {
                            if (tmpB[i] == s) indexB = i;
                        }
                        mpz_class beta = indexA;
                        mpz_class y = indexB;


                        mpz_class t = beta + y * W;
                        if (!(E.doubleAndAdd(P, (p+1+t)) == O))
                        {
                            t = beta - y * W;
                        }
                        else if ((E.doubleAndAdd(P, (p + 1 + t)) == O))
                        {
                            t = t;
                        }
                        return p + 1 + sigma * t;
                    }
                }
            }
        }
    }

    std::ostream& operator << (std::ostream& os, const EllipticCurve& curve)
    {
        std::string outputA = "A: " + curve.a.get_str() + "\n";
        std::string outputB = "B: " + curve.b.get_str() + "\n";
        std::string outputP = "P: " + curve.p.get_str() + "\n";
        std::string outputString = outputA + outputB + outputP;
        std::cout << outputString;
        return std::cout;
    }
}

