#ifndef POLY_H
#define POLY_H

struct Poly
{
public:
    double c;
    double d[6];
    Poly();
    Poly(Poly const& a);
    Poly(double a);
    Poly& operator=(Poly const& a);
    Poly& operator=(double a);
    Poly& operator+=(Poly const& a);
    Poly& operator+=(double a);
    Poly operator+(Poly const& a) const;
    Poly operator+(double a) const;
    Poly operator-() const;
    Poly& operator-=(Poly const& a);
    Poly& operator-=(double a);
    Poly operator-(Poly const& a) const;
    Poly operator-(double a) const;
    Poly& operator*=(Poly const& a);
    Poly& operator*=(double a);
    Poly operator*(Poly const& a) const;
    Poly operator*(double a) const;
    Poly& operator/=(Poly const& a);
    Poly& operator/=(double a);
    Poly operator/(Poly const& a) const;
    Poly operator/(double a) const;
};

Poly operator+(double a, Poly const& b);
Poly operator-(double a, Poly const& b);
Poly operator*(double a, Poly const& b);
Poly operator/(double a, Poly const& b);

Poly acos(Poly const& a);
Poly asin(Poly const& a);
Poly atan(Poly const& a);
Poly cos(Poly const& a);
Poly cosh(Poly const& a);
Poly exp(Poly const& a);
Poly log(Poly const& a);
Poly log10(Poly const& a);
Poly pow(Poly const& a, double y);
Poly pow(Poly const& a, int n);
Poly sin(Poly const& a);
Poly sinh(Poly const& a);
Poly sqrt(Poly const& a);
Poly tan(Poly const& a);
Poly tanh(Poly const& a);

#endif // POLY_H

//friend TJet sqrt<> ( TJet<T> const& );
//friend TJet tan<>  ( TJet<T> const& );
//friend TJet tanh<> ( TJet<T> const& );
//friend TJet erfc<> ( TJet<T> const& );

//friend TJet<std::complex<double> > erf    ( const TJet<std::complex<double> >& );
//friend TJet<double > erf    ( const TJet<double >& );
//friend TJet<std::complex<double> > w ( const TJet<std::complex<double> >& );
