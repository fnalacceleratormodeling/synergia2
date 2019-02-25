#ifndef PERIOD_H
#define PERIOD_H

#include "synergia/bunch/diagnostics_apertures_loss.h"
#include "synergia/bunch/bunch.h"

void apply_longitudinal_periodicity (Bunch & bunch, double length);
void apply_zcut (Bunch & bunch, double length, Diagnostics_loss_sptr diag_loss_sptr=boost::shared_ptr< Diagnostics_loss>());
#endif /* PERIOD_H_ */
