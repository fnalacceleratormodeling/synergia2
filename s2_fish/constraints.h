#ifndef HAVE_CONSTRAINTS_H
#define HAVE_CONSTRAINTS_H
#include "macro_bunch_store.h"

void apply_longitudinal_periodicity_t(Macro_bunch_store &mbs);
void apply_longitudinal_periodicity_z(Macro_bunch_store &mbs, double length);
void apply_circular_aperture(Macro_bunch_store &mbs, double radius);

#endif // HAVE_CONSTRAINTS_H
