#ifndef SPACE_CHARGE_OPENPMD_WRITER_H
#define SPACE_CHARGE_OPENPMD_WRITER_H

#include <openPMD/openPMD.hpp>

#include "rectangular_grid_domain.h"
#include "synergia/bunch/bunch.h"

class Space_charge_openPMD_writer {
  public:
    Space_charge_openPMD_writer(std::string const& file = "openpmd.h5");

    void
    set_write_interval(int interval)
    {
        write_interval = interval;
    }

    int
    get_write_interval() const
    {
        return write_interval;
    }

    bool start_iteration();

    // void write_rho_2d(Distributed_rectangular_grid_sptr rho);
    // void write_rho_3d(Distributed_rectangular_grid_sptr rho);

    // 2d mesh, write Ex, Ey, and Ez
    // Ez always 0
    void write_E_2d(Rectangular_grid_domain E);

    // 3d mesh, write Ex, Ey, and Ez
    // dim = 0: En is Ex
    // dim = 1: En is Ey
    // dim = 2: En is Ez
    void write_En(Rectangular_grid_domain En, int dim);

    // position, momentum, id
    void write_particles(Bunch const& bunch);

  private:
    openPMD::Series series;

    int iteration;
    int write_interval;
    bool write;
};

#endif
