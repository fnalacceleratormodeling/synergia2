#ifndef HAVE_S2_DIAGNOSTICS_H
#define HAVE_S2_DIAGNOSTICS_H

#include "macro_bunch_store.h"
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <numpy/arrayobject.h>

void
get_spatial_means(Macro_bunch_store& mbs,
                       numeric::array& means);


void
get_spatial_means_stds(Macro_bunch_store& mbs,
                       numeric::array& means,
                       numeric::array& stds);

void
get_spatial_minmax(Macro_bunch_store& mbs,
                       numeric::array& bmin,
                       numeric::array& bmax);

void
get_moments_corrs(Macro_bunch_store& mbs,
                  numeric::array& units,
                  numeric::array& means,
                  numeric::array& mom2s,
                  numeric::array& corrs,
                  numeric::array& diagmom4s);
#endif // HAVE_S2_DIAGNOSTICS_H
