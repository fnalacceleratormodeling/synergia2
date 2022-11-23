#ifndef SYNERGIA_CONFIG_H_
#define SYNERGIA_CONFIG_H_

#ifdef SYNERGIA_HAVE_OPENPMD
#include <openPMD/openPMD.hpp>
#else
#include "synergia/utils/hdf5_file.h"
#endif

#ifdef SYNERGIA_HAVE_OPENPMD
using io_device = openPMD::Series;
#else
using io_device = Hdf5_file;
#endif

#endif
