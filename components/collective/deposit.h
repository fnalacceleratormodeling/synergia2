#ifndef DEPOSIT_H_
#define DEPOSIT_H_
#include "components/collective/rectangular_grid.h"
#include "components/bunch/bunch.h"

void
deposit_charge_rectangular(Rectangular_grid & rho, Bunch & bunch,
        bool zero_first = true);

#endif /* DEPOSIT_H_ */
