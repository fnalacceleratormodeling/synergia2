#include "space_charge_3d_open_hockney.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
#include "deposit.h"

void
Space_charge_3d_open_hockney::setup_nondoubled_communication()
{
    std::vector<int > ranks1; // ranks with data from the undoubled domain
    in_group1 = false;
    for (int rank = 0; rank < comm2.get_size(); ++rank) {
        int uppers2 = distributed_fft3d_sptr->get_uppers()[rank];
        int uppers1 = std::min(uppers2, grid_shape[0]);
        int lower = 0;
        if (uppers1 <= grid_shape[0]) {
            ranks1.push_back(rank);
            if (rank == comm2.get_rank()) {
                in_group1 = true;
            }
            int length0;
            if (rank > 0) {
                length0 = uppers1 - distributed_fft3d_sptr->get_uppers()[rank
                        - 1];
            } else {
                length0 = uppers1;
            }
            lowers1.push_back(lower);
            int total_length = length0 * grid_shape[1] * grid_shape[2];
            lengths1.push_back(total_length);
            lower += total_length;
        }
    }
    int error;
    error = MPI_Comm_group(comm2.get(), &group2);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Comm_group)");
    }
    error = MPI_Group_incl(group2, ranks1.size(), &ranks1[0], &group1);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Group_incl)");
    }
    error = MPI_Comm_create(comm2.get(), group1, &mpi_comm1);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Comm_create)");
    }
    comm1.set(mpi_comm1);
}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(
        std::vector<int > const & grid_shape, bool periodic_z,
        Commxx const& comm, double z_period, double n_sigma) :
    Collective_operator("space charge"), grid_shape(3), doubled_grid_shape(3),
            padded_grid_shape(3), periodic_z(periodic_z), comm2(comm),
            z_period(z_period), n_sigma(n_sigma), domain_fixed(false)
{
    this->grid_shape[0] = grid_shape[2];
    this->grid_shape[1] = grid_shape[1];
    this->grid_shape[2] = grid_shape[0];
    for (int i = 0; i < 3; ++i) {
        doubled_grid_shape[i] = 2 * this->grid_shape[i];
    }
    distributed_fft3d_sptr = Distributed_fft3d_sptr(new Distributed_fft3d(
            doubled_grid_shape, comm));
    padded_grid_shape = distributed_fft3d_sptr->get_padded_shape_real();
    setup_nondoubled_communication();
}

Space_charge_3d_open_hockney::Space_charge_3d_open_hockney(bool periodic_z,
        Distributed_fft3d_sptr distributed_fft3d_sptr, double z_period,
        double n_sigma) :
    Collective_operator("space charge"), grid_shape(3), doubled_grid_shape(3),
            padded_grid_shape(3), periodic_z(periodic_z),
            distributed_fft3d_sptr(distributed_fft3d_sptr), comm2(
                    distributed_fft3d_sptr->get_comm()), z_period(z_period),
            n_sigma(n_sigma), domain_fixed(false)
{
    doubled_grid_shape = distributed_fft3d_sptr->get_shape();
    for (int i = 0; i < 3; ++i) {
        grid_shape[i] = doubled_grid_shape[i] / 2;
    }
    padded_grid_shape = distributed_fft3d_sptr->get_padded_shape_real();
    setup_nondoubled_communication();
}

double
Space_charge_3d_open_hockney::get_n_sigma() const
{
    return n_sigma;
}

void
Space_charge_3d_open_hockney::set_doubled_domain()
{
    std::vector<double > doubled_size(3);
    for (int i = 0; i < 3; ++i) {
        doubled_size[i] = 2 * domain_sptr->get_physical_size()[i];
    }
    doubled_domain_sptr = Rectangular_grid_domain_sptr(
            new Rectangular_grid_domain(doubled_size,
                    domain_sptr->get_physical_offset(), doubled_grid_shape,
                    periodic_z));
}

void
Space_charge_3d_open_hockney::set_fixed_domain(
        Rectangular_grid_domain_sptr domain_sptr)
{
    if ((domain_sptr->get_grid_shape()[0] != grid_shape[0])
            || (domain_sptr->get_grid_shape()[1] != grid_shape[1])
            || (domain_sptr->get_grid_shape()[2] != grid_shape[2])) {
        throw runtime_error(
                "Space_charge_3d_open_hockney::set_fixed_domain requires a shape\nequal to that of the parent object, but with zyx ordering.");
    }
    this->domain_sptr = domain_sptr;
    set_doubled_domain();
    domain_fixed = true;
}

void
Space_charge_3d_open_hockney::update_domain(Bunch const& bunch)
{
    if (!domain_fixed) {
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
        set_doubled_domain();
    }
}

Rectangular_grid_domain_sptr
Space_charge_3d_open_hockney::get_domain_sptr() const
{
    return domain_sptr;
}

Rectangular_grid_domain_sptr
Space_charge_3d_open_hockney::get_doubled_domain_sptr() const
{
    return doubled_domain_sptr;
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
    update_domain(bunch);
    Rectangular_grid_sptr local_rho_sptr(new Rectangular_grid(domain_sptr));
    deposit_charge_rectangular_zyx(*local_rho_sptr, bunch);
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
    Rectangular_grid semiglobal_rho(local_charge_density.get_domain_sptr());
    int upper = distributed_fft3d_sptr->get_upper();
    std::vector<int > uppers(distributed_fft3d_sptr->get_uppers());
    std::vector<int > real_lengths(distributed_fft3d_sptr->get_lengths());
    if (uppers[0] > grid_shape[0]) {
        real_lengths[0] = grid_shape[0] * grid_shape[1] * grid_shape[2];
    } else {
        real_lengths[0] = uppers[0] * grid_shape[1] * grid_shape[2];
    }
    for (int i = 1; i < comm2.get_size(); ++i) {
        if (uppers[i - 1] >= grid_shape[0]) {
            real_lengths[i] == 0;
        } else {
            real_lengths[i] == (uppers[i] - uppers[i - 1]) * grid_shape[1]
                    * grid_shape[2];
        }
    }
    int real_lower;
    if (comm2.get_rank() > 0) {
        real_lower = uppers[comm2.get_rank() - 1];
    } else {
        real_lower = 0;
    }
    const double * source = local_charge_density.get_grid_points().origin();

    // The destination for Reduce_scatter is the appropriate slice of semiglobal_rho
    double * dest;
    MArray3d_ref dest_points(semiglobal_rho.get_grid_points());
    dest = dest_points.origin() + (dest_points.index_bases()[0] + real_lower)
            * dest_points.strides()[0];
    MPI_Reduce_scatter((void *) source, (void *) dest, &real_lengths[0],
            MPI_DOUBLE, MPI_SUM, comm2.get());
    Distributed_rectangular_grid_sptr rho2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, real_lower,
                    upper, distributed_fft3d_sptr->get_padded_shape_real()));
    for (int i = real_lower; i < upper; ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            for (int k = 0; k < doubled_grid_shape[2]; ++k) {
                rho2->get_grid_points()[i][j][k] = 0.0;
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
Space_charge_3d_open_hockney::get_green_fn2_pointlike()
{
    if (doubled_domain_sptr == NULL) {
        throw runtime_error(
                "Space_charge_3d_open_hockney::get_green_fn2 called before domain specified");
    }
    int lower = distributed_fft3d_sptr->get_lower();
    int upper = distributed_fft3d_sptr->get_upper();
    Distributed_rectangular_grid_sptr G2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, lower, upper,
                    distributed_fft3d_sptr->get_padded_shape_real()));

    double hx = domain_sptr->get_cell_size()[2];
    double hy = domain_sptr->get_cell_size()[1];
    double hz = domain_sptr->get_cell_size()[0];

    double G000 = (1.0 / 4.0 * pi) * (3.0 / (2.0 * (sqrt(3.0)) * sqrt(hx * hx
            + hy * hy + hz * hz)));

    const int num_images = 8;
    int mix, miy; // mirror indices for x- and y-planes
    double x, y, z, G;
    const double epsz = 1.0e-12 * hz;

    for (int iz = lower; iz < upper; ++iz) {
        if (iz > grid_shape[0]) {
            z = (doubled_grid_shape[0] - iz) * hz;
        } else {
            z = iz * hz;
        }
        for (int iy = 0; iy <= grid_shape[1]; ++iy) {
            y = iy * hy;
            miy = doubled_grid_shape[1] - iy;
            if (miy == grid_shape[1]) {
                miy = doubled_grid_shape[1]; // will get thrown away
            }
            for (int ix = 0; ix <= grid_shape[2]; ++ix) {
                x = ix * hx;
                mix = doubled_grid_shape[2] - ix;
                if (mix == grid_shape[2]) {
                    mix = doubled_grid_shape[2]; // will get thrown away
                }

                if ((x == 0.0) && (y == 0.0) && (z == 0.0)) {
                    G = G000;
                } else {
                    G = 1.0 / (4.0 * pi * sqrt(x * x + y * y + z * z));
                }
                if (periodic_z) {
                    throw std::runtime_error(
                            "Space_charge_3d_open_hockney::get_green_fn2: periodic_z not yet implemented");
                    for (int image = -num_images; image < num_images; ++image) {
                        // fill me in!
                    }
                }

                G2->get_grid_points()[iz][iy][ix] = G;
                // three mirror images
                if (miy < doubled_grid_shape[1]) {
                    G2->get_grid_points()[iz][miy][ix] = G;
                    if (mix < doubled_grid_shape[2]) {
                        G2->get_grid_points()[iz][miy][mix] = G;
                    }
                }
                if (mix < doubled_grid_shape[2]) {
                    G2->get_grid_points()[iz][iy][mix] = G;
                }
            }
        }
    }

    G2->set_normalization(1.0);

    return G2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_green_fn2_linear()
{
    if (doubled_domain_sptr == NULL) {
        throw runtime_error(
                "Space_charge_3d_open_hockney::get_green_fn2 called before domain specified");
    }
    int lower = distributed_fft3d_sptr->get_lower();
    int upper = distributed_fft3d_sptr->get_upper();
    Distributed_rectangular_grid_sptr G2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, lower, upper,
                    distributed_fft3d_sptr->get_padded_shape_real()));

    double hx = domain_sptr->get_cell_size()[2];
    double hy = domain_sptr->get_cell_size()[1];
    double hz = domain_sptr->get_cell_size()[0];

    double rr = hx * hx + hy * hy;
    double r1 = sqrt(hx * hx + hy * hy + hz * hz);
    double G000 = (2.0 / rr) * (hz * r1 + rr * log((hz + r1) / sqrt(rr)) - hz
            * hz);// average value of outer cylinder.


    const int num_images = 8;
    int mix, miy; // mirror indices for x- and y-planes
    double x, y, z, G;
    const double epsz = 1.0e-12 * hz;

    for (int iz = lower; iz < upper; ++iz) {
        if (iz > grid_shape[0]) {
            z = (doubled_grid_shape[0] - iz) * hz;
        } else {
            z = iz * hz;
        }
        for (int iy = 0; iy <= grid_shape[1]; ++iy) {
            y = iy * hy;
            miy = doubled_grid_shape[1] - iy;
            if (miy == grid_shape[1]) {
                miy = doubled_grid_shape[1]; // will get thrown away
            }
            for (int ix = 0; ix <= grid_shape[2]; ++ix) {
                x = ix * hx;
                rr = x * x + y * y;
                mix = doubled_grid_shape[2] - ix;
                if (mix == grid_shape[2]) {
                    mix = doubled_grid_shape[2]; // will get thrown away
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
                    if (mix < doubled_grid_shape[2]) {
                        G2->get_grid_points()[iz][miy][mix] = G;
                    }
                }
                if (mix < doubled_grid_shape[2]) {
                    G2->get_grid_points()[iz][iy][mix] = G;
                }
            }
        }
    }

    double normalization = 1.0 / (4.0 * pi);
    G2->set_normalization(normalization);

    return G2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_scalar_field2(
        Distributed_rectangular_grid & charge_density2,
        Distributed_rectangular_grid & green_fn2)
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
    distributed_fft3d_sptr->transform(charge_density2.get_grid_points(),
            rho2hat);
    distributed_fft3d_sptr->transform(green_fn2.get_grid_points(), G2hat);

    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < cshape[1]; ++j) {
            for (int k = 0; k < cshape[2]; ++k) {
                phi2hat[i][j][k] = rho2hat[i][j][k] * G2hat[i][j][k];
            }
        }
    }

    Distributed_rectangular_grid_sptr phi2(new Distributed_rectangular_grid(
            doubled_domain_sptr, lower, upper,
            distributed_fft3d_sptr->get_padded_shape_real()));
    distributed_fft3d_sptr->inv_transform(phi2hat, phi2->get_grid_points());

    double normalization = 1.0;
    normalization *= charge_density2.get_normalization();
    normalization *= green_fn2.get_normalization();
    std::vector<int > shape(distributed_fft3d_sptr->get_padded_shape_real());
    normalization *= 1.0 / (shape[0] * shape[1] * shape[2]);
    phi2->set_normalization(normalization);

    return phi2;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::extract_scalar_field(
        Distributed_rectangular_grid const & phi2)
{
    int lower = std::min(phi2.get_lower(), grid_shape[0]);
    int upper = std::min(phi2.get_upper(), grid_shape[0]);
    Distributed_rectangular_grid_sptr phi(new Distributed_rectangular_grid(
            domain_sptr, lower, upper));

    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < grid_shape[1]; ++j) {
            for (int k = 0; k < grid_shape[2]; ++k) {
                phi->get_grid_points()[i][j][k]
                        = phi2.get_grid_points()[i][j][k];
            }
        }
    }
    phi->set_normalization(phi2.get_normalization());
    return phi;
}

Distributed_rectangular_grid_sptr
Space_charge_3d_open_hockney::get_electric_field_component(
        Distributed_rectangular_grid const& phi, int component)
{
    Distributed_rectangular_grid_sptr En(new Distributed_rectangular_grid(
            domain_sptr, phi.get_lower(), phi.get_upper()));
    MArray3d_ref En_a(En->get_grid_points());
    int lower_limit, upper_limit;
    if (component == 0) {
        lower_limit = En->get_lower_guard();
        upper_limit = En->get_upper_guard();
    } else {
        lower_limit = 0;
        upper_limit = domain_sptr->get_grid_shape()[component];
    }
    double cell_size = domain_sptr->get_cell_size()[component];
    boost::array<MArray3d::index, 3 > center, left, right;
    for (int i = En->get_lower(); i < En->get_upper(); ++i) {
        left[0] = i;
        center[0] = i;
        right[0] = i;
        for (int j = 0; j < domain_sptr->get_grid_shape()[1]; ++j) {
            left[1] = j;
            center[1] = j;
            right[1] = j;
            for (int k = 0; k < domain_sptr->get_grid_shape()[2]; ++k) {
                left[2] = k;
                center[2] = k;
                right[2] = k;

                double delta;
                if (center[component] == lower_limit) {
                    right[component] = center[component] + 1;
                    delta = cell_size;
                } else if (center[component] == upper_limit - 1) {
                    left[component] = center[component] - 1;
                    delta = cell_size;
                } else {
                    right[component] = center[component] + 1;
                    left[component] = center[component] - 1;
                    delta = 2.0 * cell_size;
                }
                En_a(center) = (En_a(right) - En_a(left)) / delta;
            }
        }
    }
    return En;
}

Rectangular_grid_sptr
Space_charge_3d_open_hockney::get_global_electric_field_component(
        Distributed_rectangular_grid const& dist_field)
{
    Rectangular_grid_sptr global_field(new Rectangular_grid(domain_sptr));
    const int root = 0;
    int error;
    if (in_group1) {
        error = MPI_Gatherv((void *) (dist_field.get_grid_points().origin()
                + lowers1[comm1.get_rank()]), lengths1[comm1.get_rank()],
                MPI_DOUBLE, (void*) global_field->get_grid_points().origin(),
                &lengths1[0], &lowers1[0], MPI_DOUBLE, root, comm1.get());
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                    "MPI error in Space_charge_3d_open_hockney(MPI_Gatherv)");
        }

    }
    int total_length = grid_shape[0] * grid_shape[1] * grid_shape[2];
    error = MPI_Bcast(global_field->get_grid_points().origin(), total_length,
            MPI_DOUBLE, root, comm2.get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Bcast)");
    }
    return global_field;
}

inline double
interpolate_rectangular_zyx(double x, double y, double z,
        Rectangular_grid const& f)
{
    // tri-linear interpolation
    int ix, iy, iz;
    double offx, offy, offz;
    f.get_domain_sptr()->get_leftmost_indices_offsets(z, y, x, iz, iy, ix,
            offz, offy, offx);
    MArray3d_ref a(f.get_grid_points());
    double val;
    if ((iz < 0) || (iz >= f.get_domain_sptr()->get_grid_shape()[0] - 1) || (iy
            < 0) || (iy >= f.get_domain_sptr()->get_grid_shape()[1] - 1) || (ix
            < 0) || (ix >= f.get_domain_sptr()->get_grid_shape()[2] - 1)) {
        val = 0.0;
    } else {
        val = ((1.0 - offz) * (1.0 - offy) * (1.0 - offx) * a[iz][iy][ix]
                + offz * (1.0 - offy) * (1.0 - offx) * a[iz + 1][iy][ix] + (1.0
                - offz) * offy * (1.0 - offx) * a[iz][iy + 1][ix]
                + (1.0 - offz) * (1.0 - offy) * offx * a[iz][iy][ix + 1] + offz
                * offy * (1.0 - offx) * a[iz + 1][iy + 1][ix] + offz * (1.0
                - offy) * offx * a[iz + 1][iy][ix + 1] + (1.0 - offz) * offy
                * offx * a[iz][iy + 1][ix + 1] + offz * offy * offx
                * a[iz + 1][iy + 1][ix + 1]);
    }
    return val;
}

void
Space_charge_3d_open_hockney::apply_kick(Bunch & bunch,
        Rectangular_grid const& En, double delta_tau, int component)
{
    // $\delta \vec{p} = \vec{F} \delta t = q \delta t \vec{E}$
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    double factor = q * delta_tau * En.get_normalization();
    int ps_component = 2 * component + 1;
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double x = bunch.get_local_particles()[part][Bunch::x];
        double y = bunch.get_local_particles()[part][Bunch::y];
        double z = bunch.get_local_particles()[part][Bunch::z];
        double field_val = interpolate_rectangular_zyx(x, y, z, En);
        bunch.get_local_particles()[part][ps_component] += factor * field_val;
    }
}

void
Space_charge_3d_open_hockney::apply(Bunch & bunch, double time_step,
        Step & step)
{
    Rectangular_grid_sptr local_rho(get_local_charge_density(bunch)); // [C/m^3]
    Distributed_rectangular_grid_sptr rho2(get_global_charge_density2(
            *local_rho)); // [C/m^3]
    local_rho.reset();
    Distributed_rectangular_grid_sptr G2(get_green_fn2_pointlike()); // [1/m^3]
    Distributed_rectangular_grid_sptr phi2(get_scalar_field2(*rho2, *G2)); // [V]
    rho2.reset();
    G2.reset();
    Distributed_rectangular_grid_sptr phi(extract_scalar_field(*phi2));
    phi2.reset();
    phi->fill_guards(comm1);
    for (int component = 0; component < 3; ++component) {
        Distributed_rectangular_grid_sptr local_En(
                get_electric_field_component(*phi, component)); // [V/m]
        Rectangular_grid_sptr
                En(get_global_electric_field_component(*local_En)); // [V/m]
        apply_kick(bunch, *En, time_step, component);
    }
}

Space_charge_3d_open_hockney::~Space_charge_3d_open_hockney()
{
    int error;
    error = MPI_Comm_free(&mpi_comm1);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Comm_free)");
    }
    error = MPI_Group_free(&group1);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Group_free(1))");
    }
    error = MPI_Group_free(&group2);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_3d_open_hockney(MPI_Group_free(2))");
    }

}

