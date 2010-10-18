#ifndef DEPOSIT_H_
#define DEPOSIT_H_
#include "synergia/collective/rectangular_grid.h"
#include "synergia/bunch/bunch.h"

void
deposit_charge_rectangular(Rectangular_grid & rho, Bunch const& bunch,
        bool zero_first = true);

#endif /* DEPOSIT_H_ */
