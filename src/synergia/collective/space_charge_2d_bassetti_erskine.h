#ifndef SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_
#define SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_

#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"

class Space_charge_2d_bassetti_erskine : public Collective_operator
{
public:
    Space_charge_2d_bassetti_erskine();
    virtual void
    apply(Bunch & bunch, double time_step, Step & step);
    virtual
    ~Space_charge_2d_bassetti_erskine();
};

typedef boost::shared_ptr<Space_charge_2d_bassetti_erskine >
        Space_charge_2d_bassetti_erskine_sptr;

#endif /* SPACE_CHARGE_2D_BASSETTI_ERSKINE_H_ */
