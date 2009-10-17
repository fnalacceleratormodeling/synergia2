#include "diagnostics.h"

//class Diagnostics
//{
//protected:
//    double s;
//    MArray1d mean;
//    MArray1d std;

void
Diagnostics::update_mean(Bunch const& bunch)
{
}

void
Diagnostics::update_std(Bunch const& bunch)
{
}

Diagnostics::Diagnostics() :
    mean(boost::extents[6]), std(boost::extents[6])
{
}

Diagnostics::Diagnostics(Bunch const& bunch, double s):
    mean(boost::extents[6]), std(boost::extents[6])
{
    update(bunch,s);
}

void
Diagnostics::update(Bunch const& bunch, double s)
{
    this->s = s;
    update_mean(bunch);
    update_std(bunch);
}

double
Diagnostics::get_s() const
{
    return s;
}

Const_MArray1d_ref
Diagnostics::get_mean() const
{
}

Const_MArray1d_ref
Diagnostics::get_std() const
{
}

Diagnostics::~Diagnostics()
{
}

Diagnostics_full2::Diagnostics_full2()
{
}

Diagnostics_full2::Diagnostics_full2(Bunch const& bunch, double s)
{
}

void
Diagnostics_full2::update(Bunch const& bunch, double s)
{
}

Const_MArray2d_ref
Diagnostics_full2::get_mom2() const
{
}

Const_MArray2d_ref
Diagnostics_full2::get_corr() const
{
}

double
Diagnostics_full2::get_emitx() const
{
}

double
Diagnostics_full2::get_emity() const
{
}

double
Diagnostics_full2::get_emitz() const
{
}

double
Diagnostics_full2::get_emitxyz() const
{
}

Diagnostics_full2::~Diagnostics_full2()
{
}
