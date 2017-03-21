#include "poly.h"
#include <cmath>

Poly::Poly() :
    c(0.0),
    d()
{
    d[0] = 0.0;
    d[1] = 0.0;
    d[2] = 0.0;
    d[3] = 0.0;
    d[4] = 0.0;
    d[5] = 0.0;
}

Poly::Poly(Poly const& a) :
    c(a.c),
    d()
{
    d[0] = a.d[0];
    d[1] = a.d[1];
    d[2] = a.d[2];
    d[3] = a.d[3];
    d[4] = a.d[4];
    d[5] = a.d[5];
}

Poly::Poly(double a) :
    c(a),
    d()
{
    d[0] = 0.0;
    d[1] = 0.0;
    d[2] = 0.0;
    d[3] = 0.0;
    d[4] = 0.0;
    d[5] = 0.0;
}

Poly& Poly::operator=(Poly const& a)
{
    c = a.c;
    d[0] = a.d[0];
    d[1] = a.d[1];
    d[2] = a.d[2];
    d[3] = a.d[3];
    d[4] = a.d[4];
    d[5] = a.d[5];

    return *this;
}

Poly& Poly::operator=(double a)
{
    c = a;
    d[0] = 0.0;
    d[1] = 0.0;
    d[2] = 0.0;
    d[3] = 0.0;
    d[4] = 0.0;
    d[5] = 0.0;

    return *this;
}

Poly& Poly::operator+=(Poly const& a)
{
    c += a.c;
    d[0] += a.d[0];
    d[1] += a.d[1];
    d[2] += a.d[2];
    d[3] += a.d[3];
    d[4] += a.d[4];
    d[5] += a.d[5];

    return *this;
}

Poly& Poly::operator+=(double a)
{
    c += a;

    return *this;
}

Poly Poly::operator+(Poly const& a) const
{
    Poly retval(*this);
    retval += a;

    return retval;
}

Poly Poly::operator+(double a) const
{
    Poly retval(*this);
    retval += a;

    return retval;
}

Poly Poly::operator-() const
{
    Poly retval(*this);

    retval.c = -c;
    retval.d[0] = -d[0];
    retval.d[1] = -d[1];
    retval.d[2] = -d[2];
    retval.d[3] = -d[3];
    retval.d[4] = -d[4];
    retval.d[5] = -d[5];

    return retval;
}

Poly& Poly::operator-=(Poly const& a)
{
    c -= a.c;
    d[0] -= a.d[0];
    d[1] -= a.d[1];
    d[2] -= a.d[2];
    d[3] -= a.d[3];
    d[4] -= a.d[4];
    d[5] -= a.d[5];

    return *this;
}

Poly& Poly::operator-=(double a)
{
    c -= a;

    return *this;
}

Poly Poly::operator-(Poly const& a) const
{
    Poly retval(*this);
    retval -= a;

    return retval;
}

Poly Poly::operator-(double a) const
{
    Poly retval(*this);
    retval -= a;

    return retval;
}

Poly& Poly::operator*=(Poly const& a)
{
    d[0] = c * a.d[0] + d[0] * a.c;
    d[1] = c * a.d[1] + d[1] * a.c;
    d[2] = c * a.d[2] + d[2] * a.c;
    d[3] = c * a.d[3] + d[3] * a.c;
    d[4] = c * a.d[4] + d[4] * a.c;
    d[5] = c * a.d[5] + d[5] * a.c;
    c *= a.c;

    return *this;
}

Poly& Poly::operator*=(double a)
{
    c *= a;
    d[0] *= a;
    d[1] *= a;
    d[2] *= a;
    d[3] *= a;
    d[4] *= a;
    d[5] *= a;

    return *this;
}

Poly Poly::operator*(Poly const& a) const
{
    Poly retval(*this);
    retval *= a;

    return retval;
}

Poly Poly::operator*(double a) const
{
    Poly retval(*this);
    retval *= a;

    return retval;
}

Poly& Poly::operator/=(Poly const& a)
{
    const double ac2(a.c * a.c);
    d[0] = (- c * a.d[0] + d[0] * a.c)/ac2;
    d[1] = (- c * a.d[1] + d[1] * a.c)/ac2;
    d[2] = (- c * a.d[2] + d[2] * a.c)/ac2;
    d[3] = (- c * a.d[3] + d[3] * a.c)/ac2;
    d[4] = (- c * a.d[4] + d[4] * a.c)/ac2;
    d[5] = (- c * a.d[5] + d[5] * a.c)/ac2;
    c /= a.c;

    return *this;
}

Poly& Poly::operator/=(double a)
{
    c /= a;
    d[0] /= a;
    d[1] /= a;
    d[2] /= a;
    d[3] /= a;
    d[4] /= a;
    d[5] /= a;

    return *this;
}

Poly Poly::operator/(Poly const& a) const
{
    Poly retval(*this);
    retval /= a;

    return retval;
}

Poly Poly::operator/(double a) const
{
    Poly retval(*this);
    retval /= a;

    return retval;
}

Poly operator+(double a, Poly const& b)
{
    return b + a;
}

Poly operator-(double a, Poly const& b)
{
    return -(b - a);
}

Poly operator*(double a, Poly const& b)
{
    return b * a;
}

Poly operator/(double a, Poly const& b)
{
    Poly polya(a);

    return polya/b;
}

Poly acos(Poly const& a)
{
    Poly retval;
    retval.c = std::acos(a.c);
    double dfdx = -1.0/std::sqrt((1.0 - a.c)*(1.0 + a.c));
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly asin(Poly const& a)
{
    Poly retval;
    retval.c = std::asin(a.c);
    double dfdx = 1.0/std::sqrt((1.0 - a.c)*(1.0 + a.c));
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly atan(Poly const& a)
{
    Poly retval;
    retval.c = std::atan(a.c);
    double dfdx = 1.0/(1.0 + a.c * a.c);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly cos(Poly const& a)
{
    Poly retval;
    retval.c = std::cos(a.c);
    double dfdx = -std::sin(a.c);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly cosh(Poly const& a)
{
    Poly retval;
    retval.c = std::cosh(a.c);
    double dfdx = std::sinh(a.c);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly exp(Poly const& a)
{
    Poly retval;
    retval.c = std::exp(a.c);
    double dfdx = retval.c;
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly log(Poly const& a)
{
    Poly retval;
    retval.c = std::log(a.c);
    double dfdx = 1.0/a.c;
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly log10(Poly const& a)
{
    Poly retval;
    retval.c = std::log10(a.c);
    double dfdx = 1.0/(std::log(10.0)*a.c);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly pow(Poly const& a, double y)
{
    Poly retval;
    retval.c = std::pow(a.c, y);
    double dfdx = y * std::pow(a.c, y - 1.0);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly pow(Poly const& a, int n)
{
    Poly retval;
    retval.c = std::pow(a.c, n);
    double dfdx = n * std::pow(a.c, n - 1);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly sin(Poly const& a)
{
    Poly retval;
    retval.c = std::sin(a.c);
    double dfdx = std::cos(a.c);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly sinh(Poly const& a)
{
    Poly retval;
    retval.c = std::sinh(a.c);
    double dfdx = std::cosh(a.c);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly sqrt(Poly const& a)
{
    Poly retval;
    retval.c = std::sqrt(a.c);
    double dfdx = 0.5/std::sqrt(a.c);
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly tan(Poly const& a)
{
    Poly retval;
    retval.c = std::tan(a.c);
    double dfdx = 1.0/(std::cos(a.c)*std::cos(a.c));
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}

Poly tanh(Poly const& a)
{
    Poly retval;
    retval.c = std::tanh(a.c);
    double dfdx = 1.0/(std::cosh(a.c)*std::cosh(a.c));;
    for(unsigned int i = 0; i < 6; ++i) {
        retval.d[i] = dfdx * a.d[i];
    }

    return retval;
}
