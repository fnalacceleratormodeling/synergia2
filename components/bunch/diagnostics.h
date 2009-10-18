#ifndef DIAGNOSTICS_H_
#define DIAGNOSTICS_H_

#include "components/bunch/bunch.h"

class Diagnostics
{
protected:
    double s;
    MArray1d mean;
    MArray1d std;
    virtual void
    update_mean(Bunch const& bunch);
    virtual void
    update_std(Bunch const& bunch);
public:
    Diagnostics();
    Diagnostics(Bunch const& bunch, double s);
    virtual void
    update(Bunch const& bunch, double s);
    virtual double
    get_s() const;
    virtual Const_MArray1d_ref
    get_mean() const;
    virtual Const_MArray1d_ref
    get_std() const;
    virtual
    ~Diagnostics();
};

class Diagnostics_full2 : public Diagnostics
{
protected:
    MArray2d mom2;
    MArray2d corr;
    double emitx, emity, emitz, emitxy, emitxyz;
    virtual void
    update_full2(Bunch const& bunch);
    virtual void
    update_emittances();
public:
    Diagnostics_full2();
    Diagnostics_full2(Bunch const& bunch, double s);
    virtual void
    update(Bunch const& bunch, double s);
    virtual Const_MArray2d_ref
    get_mom2() const;
    virtual Const_MArray2d_ref
    get_corr() const;
    virtual double
    get_emitx() const;
    virtual double
    get_emity() const;
    virtual double
    get_emitz() const;
    virtual double
    get_emitxy() const;
    virtual double
    get_emitxyz() const;
    virtual
    ~Diagnostics_full2();
};

#endif /* DIAGNOSTICS_H_ */
