#include "space_charge_3d_open_hockney.h"
#include "components/bunch/diagnostics.h"
#include "components/foundation/math_constants.h"
#include "deposit.h"

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
        std::vector<int > const & grid_shape, bool periodic_z,
        Commxx const& comm, double z_period, double n_sigma) :
    Collective_operator("space charge"), grid_shape(grid_shape),
            doubled_grid_shape(3), periodic_z(periodic_z), comm(comm),
            z_period(z_period), n_sigma(n_sigma)
{
    for (int i = 0; i < 3; ++i) {
        doubled_grid_shape[i] = 2 * grid_shape[i];
    }
    distributed_fft3d_sptr = Distributed_fft3d_sptr(new Distributed_fft3d(
            doubled_grid_shape, comm));
}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(bool periodic_z,
        Distributed_fft3d_sptr const& distributed_fft3d_sptr, double z_period,
        double n_sigma) :
    Collective_operator("space charge"), grid_shape(3), periodic_z(periodic_z),
            distributed_fft3d_sptr(distributed_fft3d_sptr), comm(
                    distributed_fft3d_sptr->get_comm()), z_period(z_period),
            n_sigma(n_sigma)
{
    grid_shape = distributed_fft3d_sptr->get_shape();
}

void
Space_charge_3d_open_hockney::update_domain(Bunch const& bunch)
{
    Diagnostics diagnostics(bunch);
    std::vector<double > size(3);
    std::vector<double > offset(3);
    offset[0] = diagnostics.get_mean()[Bunch::z];
    size[0] = n_sigma * diagnostics.get_std()[Bunch::z];
    offset[1] = diagnostics.get_mean()[Bunch::y];
    size[1] = n_sigma * diagnostics.get_std()[Bunch::y];
    offset[2] = diagnostics.get_mean()[Bunch::x];
    size[2] = n_sigma * diagnostics.get_std()[Bunch::x];
    domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
            size, offset, grid_shape, periodic_z));
    std::vector<double > doubled_size(3);
    for (int i = 0; i < 3; ++i) {
        doubled_size[i] = 2 * size[i];
    }
    doubled_domain_sptr = Rectangular_grid_domain_sptr(
            new Rectangular_grid_domain(doubled_size, offset,
                    doubled_grid_shape, periodic_z));
}

Rectangular_grid_domain_sptr
Space_charge_3d_open_hockney::get_domain_sptr()
{
    return domain_sptr;
}

Rectangular_grid_domain_sptr
Space_charge_3d_open_hockney::get_doubled_domain_sptr()
{
    return doubled_domain_sptr;
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
    update_domain(bunch);
    Rectangular_grid_sptr local_rho_sptr(new Rectangular_grid(domain_sptr));
    deposit_charge_rectangular(*local_rho_sptr, bunch);
    return local_rho_sptr;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_charge_density2(
        Rectangular_grid const& local_charge_density)
{
    // jfa: here is where we do something complicated, but (potentially) efficient
    // in calculating a version of the charge density that is just global enough
    // to fill in the doubled global charge density

    // semiglobal_rho stores the portion of global charge density needed on each
    // processor. It has to have the same shape as the charge density in order
    // to work with MPI_Reduce_scatter.
    Rectangular_grid semiglobal_rho(
            local_charge_density.get_domain_sptr());
    int upper = distributed_fft3d_sptr->get_upper();
    std::vector<int > uppers(distributed_fft3d_sptr->get_uppers());
    std::vector<int > real_lengths(distributed_fft3d_sptr->get_lengths());
    for (int i = 1; i < comm.get_size(); ++i) {
        if (uppers[i - 1] >= grid_shape[0]) {
            real_lengths[i] == 0;
        }
    }
    int real_lower = upper - real_lengths[comm.get_rank()];
    void * source;
    source = (void *) local_charge_density.get_grid_points().origin();

    // The destination for Reduce_scatter is the appropriate slice of semiglobal_rho
    void * dest;
    MArray3d_ref dest_points(semiglobal_rho.get_grid_points());
    dest = dest_points.origin() + (dest_points.index_bases()[0] + real_lower)
            * dest_points.strides()[0];
    MPI_Reduce_scatter(source, dest, &real_lengths[0], MPI_DOUBLE, MPI_SUM,
            comm.get());

    Distributed_rectangular_grid_sptr rho2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, real_lower,
                    upper));

    for (int i = real_lower; i < upper; ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            for (int k = 0; k < doubled_grid_shape[2]; ++k) {
                semiglobal_rho.get_grid_points()[i][j][k] = 0.0;
            }
        }
    }
    int nonzero_i_max = (upper < grid_shape[0]) ? upper : grid_shape[0];
    for (int i = real_lower; i < nonzero_i_max; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                rho2->get_grid_points()[i][j][k]
                        = semiglobal_rho.get_grid_points()[i][j][k];
            }
        }
    }

    return rho2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_green_fn2()
{
    int lower = distributed_fft3d_sptr->get_upper();
    int upper = distributed_fft3d_sptr->get_upper();
    Distributed_rectangular_grid_sptr G2 =
            Distributed_rectangular_grid_sptr(new Distributed_rectangular_grid(
                    doubled_domain_sptr, lower, upper));

    double hx = domain_sptr->get_cell_size()[0];
    double hy = domain_sptr->get_cell_size()[1];
    double hz = domain_sptr->get_cell_size()[2];

    double rr = hx * hx + hy * hy;
    double r1 = sqrt(hx * hx + hy * hy + hz * hz);

    const int num_images = 8;
    int mix, miy; // mirror indices for x- and y-planes
    double x, y, z, G, G000;
    const double epsz = 1.0e-12 * hz;

    for (int iz = lower; iz < upper; ++iz) {
        if (iz > grid_shape[2]) {
            z = (doubled_grid_shape[2] - iz) * hz;
        } else {
            z = iz * hz;
        }
        for (int iy = 0; iy <= grid_shape[1]; ++iy) {
            y = iy * hy;
            miy = doubled_grid_shape[1] - iy;
            if (miy == grid_shape[1]) {
                miy = doubled_grid_shape[1]; // will get thrown away
            }
            for (int ix = 0; ix <= grid_shape[0]; ++ix) {
                x = ix * hx;
                rr = x * x + y * y;
                mix = doubled_grid_shape[0] - ix;
                if (mix == grid_shape[0]) {
                    mix = doubled_grid_shape[0]; // will get thrown away
                }

                G = 2.0 * sqrt(rr + z * z) - sqrt(rr + (z - hz) * (z - hz))
                        - sqrt(rr + (z + hz) * (z + hz));
                double T1, T2, r1, r2;
                if (z < -hz) {
                    r1 = (sqrt((z - hz) * (z - hz) + rr) - z + hz) / (sqrt(z
                            * z + rr) - z);
                    T1 = (hz - z) * log(r1);
                    r2 = (sqrt(z * z + rr) - z) / (sqrt((z + hz) * (z + hz)
                            + rr) - z - hz);
                    T2 = (hz + z) * log(r2);
                    G += T1 + T2;
                } else if (fabs(z + hz) < epsz) {
                    r1 = (sqrt((z - hz) * (z - hz) + rr) - z + hz) / (sqrt(z
                            * z + rr) - z);
                    T1 = (hz - z) * log(r1);
                    G += T1;
                } else if (fabs(z) < epsz) {
                    if (fabs(x) + fabs(y) < 2. * epsz) {
                        G += hz * G000;
                    } /* T1+T2 in fact */else {
                        r1 = (sqrt(hz * hz + rr) + hz) / sqrt(rr);
                        G += 2. * hz * log(r1);
                    }
                } else if (fabs(z - hz) < epsz) {
                    r1 = (sqrt((z + hz) * (z + hz) + rr) + z + hz) / (sqrt(z
                            * z + rr) + z);
                    T1 = (hz + z) * log(r1);
                    G += T1;
                } else if (z > hz) {
                    r1 = (sqrt(z * z + rr) + z) / (sqrt((z - hz) * (z - hz)
                            + rr) + z - hz);
                    T1 = (hz - z) * log(r1);
                    r2 = (sqrt((z + hz) * (z + hz) + rr) + z + hz) / (sqrt(z
                            * z + rr) + z);
                    T2 = (hz + z) * log(r2);
                    G += T1 + T2;
                } else {
                    throw std::runtime_error(
                            "Space_charge_3d_open_hockney::get_green_fn2 error1");
                }

                if (periodic_z) {
                    throw std::runtime_error(
                            "Space_charge_3d_open_hockney::get_green_fn2: periodic_z not yet implemented");
                    for (int image = -num_images; image < num_images; ++image) {
                        if (image != 0) {
                            double z_image = z + image * z_period;

                            if (z_image < -hz) {
                                r1
                                        = (sqrt((z_image - hz) * (z_image - hz)
                                                + rr) - z_image + hz) / (sqrt(
                                                z_image * z_image + rr)
                                                - z_image);
                                T1 = (hz - z_image) * log(r1);
                                r2 = (sqrt(z_image * z_image + rr) - z_image)
                                        / (sqrt((z_image + hz) * (z_image + hz)
                                                + rr) - z_image - hz);
                                T2 = (hz + z_image) * log(r2);
                                G += T1 + T2;
                            }

                            else if (fabs(z_image + hz) < epsz) {
                                r1
                                        = (sqrt((z_image - hz) * (z_image - hz)
                                                + rr) - z_image + hz) / (sqrt(
                                                z_image * z_image + rr)
                                                - z_image);
                                T1 = (hz - z_image) * log(r1);
                                G += T1;
                            }

                            else if (fabs(z_image) < epsz) {
                                if (fabs(x) + fabs(y) < 2. * epsz) {
                                    G += hz * G000;
                                } // T1+T2 in fact
                                else {
                                    r1 = (sqrt(hz * hz + rr) + hz) / sqrt(rr);
                                    G += 2. * hz * log(r1);
                                }
                            } else if (fabs(z_image - hz) < epsz) {
                                r1
                                        = (sqrt((z_image + hz) * (z_image + hz)
                                                + rr) + z_image + hz) / (sqrt(
                                                z_image * z_image + rr)
                                                + z_image);
                                T1 = (hz + z_image) * log(r1);
                                G += T1;
                            } else if (z_image > hz) {
                                r1 = (sqrt(z_image * z_image + rr) + z_image)
                                        / (sqrt((z_image - hz) * (z_image - hz)
                                                + rr) + z_image - hz);
                                T1 = (hz - z_image) * log(r1);
                                r2
                                        = (sqrt((z_image + hz) * (z_image + hz)
                                                + rr) + z_image + hz) / (sqrt(
                                                z_image * z_image + rr)
                                                + z_image);
                                T2 = (hz + z_image) * log(r2);
                                G += T1 + T2;
                            } else {
                                throw std::runtime_error(
                                        "Space_charge_3d_open_hockney::get_green_fn2 error2");
                            }

                        }
                    }
                }

                G2->get_grid_points()[iz][iy][ix] = G;
                // three mirror images
                if (miy < doubled_grid_shape[1]) {
                    G2->get_grid_points()[iz][miy][ix] = G;
                    if (mix < doubled_grid_shape[0]) {
                        G2->get_grid_points()[iz][miy][mix] = G;
                    }
                }
                if (mix < doubled_grid_shape[0]) {
                    G2->get_grid_points()[iz][iy][mix] = G;
                }
            }
        }
    }

    double scale = 1.0 / (4.0 * constants::pi * hz * hz);
    for (int iz = lower; iz < upper; ++iz) {
        for (int iy = 0; iy < grid_shape[2]; ++iy) {
            for (int ix = 0; ix < grid_shape[0]; ++ix) {
                G2->get_grid_points()[iz][iy][ix] *= scale;
            }
        }
    }

    return G2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_scalar_field2(
        Distributed_rectangular_grid_sptr & charge_density2,
        Distributed_rectangular_grid_sptr & green_fn2)
{
    std::vector<int >
            cshape(distributed_fft3d_sptr->get_padded_shape_complex());
    int lower = distributed_fft3d_sptr->get_lower();
    int upper = distributed_fft3d_sptr->get_upper();
    MArray3dc rho2hat(
            boost::extents[extent_range(lower, upper)][cshape[1]][cshape[2]]);
    MArray3dc G2hat(
            boost::extents[extent_range(lower, upper)][cshape[1]][cshape[2]]);
    MArray3dc phi2hat(
            boost::extents[extent_range(lower, upper)][cshape[1]][cshape[2]]);

    distributed_fft3d_sptr->transform(charge_density2->get_grid_points(),
            rho2hat);
    distributed_fft3d_sptr->transform(green_fn2->get_grid_points(), G2hat);

    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < cshape[1]; ++j) {
            for (int k = 0; k < cshape[2]; ++k) {
                phi2hat[i][j][k] = rho2hat[i][j][k] * G2hat[i][j][k];
            }
        }
    }

    Distributed_rectangular_grid_sptr phi2(new Distributed_rectangular_grid(
            doubled_domain_sptr, lower, upper));
    distributed_fft3d_sptr->inv_transform(phi2hat,phi2->get_grid_points());

    return phi2;
}

void
Space_charge_3d_open_hockney::apply(Bunch & bunch, Operators & step_operators)
{
    Rectangular_grid_sptr local_rho(get_local_charge_density(bunch));
    Distributed_rectangular_grid_sptr rho2(
            get_global_charge_density2(*local_rho));
    Distributed_rectangular_grid_sptr G2(get_green_fn2());
    Distributed_rectangular_grid_sptr phi2(get_scalar_field2(rho2, G2));
}

Space_charge_3d_open_hockney::~Space_charge_3d_open_hockney()
{

}

