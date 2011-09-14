#include "space_charge_2d_open_hockney.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/foundation/math_constants.h"
using mconstants::pi;
#include "synergia/foundation/physical_constants.h"
using pconstants::epsilon0;
#include "deposit.h"
#include "interpolate_rectangular_zyx.h"
#include "synergia/utils/multi_array_offsets.h"
#include "synergia/utils/simple_timer.h"

void
Space_charge_2d_open_hockney::setup_nondoubled_communication()
{
    std::vector<int > ranks1; // ranks with data from the undoubled domain
    in_group1 = false;
    int lower = 0;
    for (int rank = 0; rank < comm2.get_size(); ++rank) {
        int uppers2 = distributed_fft2d_sptr->get_uppers()[rank];
        int uppers1 = std::min(uppers2, doubled_grid_shape[0]);
        int length0;
        if (rank > 0) {
            length0 = uppers2 - distributed_fft2d_sptr->get_uppers()[rank - 1];
        } else {
            length0 = uppers2;
        }
        if (length0 > 0) {
            ranks1.push_back(rank);
            if (rank == comm2.get_rank()) {
                in_group1 = true;
            }
            lowers1.push_back(lower);
            int total_length = length0 * doubled_grid_shape[1];
            lengths1.push_back(total_length);
            lower += total_length;
        }
    }
    int error;
    error = MPI_Comm_group(comm2.get(), &group2);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney(MPI_Comm_group)");
    }
    error = MPI_Group_incl(group2, ranks1.size(), &ranks1[0], &group1);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney(MPI_Group_incl)");
    }
    error = MPI_Comm_create(comm2.get(), group1, &mpi_comm1);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney(MPI_Comm_create)");
    }
    comm1.set(mpi_comm1);

    std::vector<int > real_uppers(distributed_fft2d_sptr->get_uppers());
    real_lengths = distributed_fft2d_sptr->get_lengths();
    real_lengths_1d = distributed_fft2d_sptr->get_lengths_1d();
    for (int i = 0; i < comm2.get_size(); ++i) {
        if (i == 0) {
            real_lengths[0] = real_uppers[0] * doubled_grid_shape[1]; 
            real_lengths_1d[0] = doubled_grid_shape[2];
        } else {
            real_lengths[i] = (real_uppers[i] - real_uppers[i - 1])
                    * doubled_grid_shape[1];
            real_lengths_1d[i] = doubled_grid_shape[2];
        }
    }
    int my_rank = comm2.get_rank();
    if (my_rank > 0) {
        real_lower = real_uppers[my_rank - 1];
    } else {
        real_lower = 0;
    }
    real_upper = real_uppers[my_rank];
    real_length = real_lengths[my_rank];
    if (my_rank > 0) {
        doubled_lower = distributed_fft2d_sptr->get_uppers()[my_rank - 1];
    } else {
        doubled_lower = 0;
    }
    doubled_upper = distributed_fft2d_sptr->get_uppers()[my_rank];
    real_doubled_lower = std::min(doubled_lower, grid_shape[0]);
    real_doubled_upper = std::min(doubled_upper, grid_shape[0]);
}

void
Space_charge_2d_open_hockney::setup_default_options()
{
    set_green_fn_type(pointlike);
    set_charge_density_comm(charge_allreduce);
    set_e_force_comm(e_force_allreduce);
}

Space_charge_2d_open_hockney::Space_charge_2d_open_hockney(Commxx const& comm,
        std::vector<int > const & grid_shape, bool periodic_z, double z_period, 
        bool grid_entire_period, double n_sigma) :
    Collective_operator("space charge 3D open hockney"), comm2(comm),
            grid_shape(3), doubled_grid_shape(3), periodic_z(periodic_z), 
            z_period(z_period), grid_entire_period(grid_entire_period), 
            n_sigma(n_sigma), domain_fixed(false), have_domains(false)
{
    if (this->periodic_z) {
        throw std::runtime_error(
                "Space_charge_2d_open_hockney: periodic_z must be FALSE");
    }
    // no xyz->zyx transformation in 2D version
    for (int i = 0; i < 3; ++i) {
        this->grid_shape[i] = grid_shape[i];
        doubled_grid_shape[i] = 2 * this->grid_shape[i];
    }
    doubled_grid_shape[2] = this->grid_shape[2]; 
    distributed_fft2d_sptr = Distributed_fft2d_sptr(new Distributed_fft2d(
            doubled_grid_shape, comm));
    setup_nondoubled_communication();
    setup_default_options();
}

Space_charge_2d_open_hockney::Space_charge_2d_open_hockney(
        Distributed_fft2d_sptr distributed_fft2d_sptr, bool periodic_z, 
        double z_period, bool grid_entire_period, double n_sigma) :
    Collective_operator("space charge"), grid_shape(3), doubled_grid_shape(3),
            comm2(distributed_fft2d_sptr->get_comm()),
            distributed_fft2d_sptr(distributed_fft2d_sptr), 
            periodic_z(periodic_z), z_period(z_period), 
            grid_entire_period(grid_entire_period), n_sigma(n_sigma), 
            domain_fixed(false), have_domains(false)
{
    if (this->periodic_z) {
        throw std::runtime_error(
                "Space_charge_2d_open_hockney: periodic_z must be FALSE");
    }
    doubled_grid_shape = distributed_fft2d_sptr->get_shape();
    for (int i = 0; i < 2; ++i) {
        grid_shape[i] = doubled_grid_shape[i] / 2;
    }
    grid_shape[2] = doubled_grid_shape[2];
    setup_nondoubled_communication();
    setup_default_options();
}

double
Space_charge_2d_open_hockney::get_n_sigma() const
{
    return n_sigma;
}

void
Space_charge_2d_open_hockney::set_green_fn_type(Green_fn_type green_fn_type)
{
    switch (green_fn_type) {
    case pointlike:
        break;
    case bruteforce:
        break;
    default:
        throw runtime_error(
                "Space_charge_2d_open_hockney::set_green_fn_type: invalid green_fn_type");
    }
    this->green_fn_type = green_fn_type;
}

Space_charge_2d_open_hockney::Green_fn_type
Space_charge_2d_open_hockney::get_green_fn_type() const
{
    return green_fn_type;
}

void
Space_charge_2d_open_hockney::set_charge_density_comm(
        Charge_density_comm charge_density_comm)
{
    switch (charge_density_comm) {
    case reduce_scatter:
        break;
    case charge_allreduce:
        break;
    default:
        throw runtime_error(
                "Space_charge_2d_open_hockney::set_charge_density_comm: invalid charge_density_comm");
    }
    this->charge_density_comm = charge_density_comm;
}

Space_charge_2d_open_hockney::Charge_density_comm
Space_charge_2d_open_hockney::get_charge_density_comm() const
{
    return charge_density_comm;
}

void
Space_charge_2d_open_hockney::set_e_force_comm(E_force_comm e_force_comm)
{
    switch (e_force_comm) {
    case gatherv_bcast:
        break;
    case allgatherv:
        break;
    case e_force_allreduce:
        break;
    default:
        throw runtime_error(
                "Space_charge_2d_open_hockney::set_e_force_comm: invalid e_force_comm");
    }

    this->e_force_comm = e_force_comm;
}

Space_charge_2d_open_hockney::E_force_comm
Space_charge_2d_open_hockney::get_e_force_comm() const
{
    return e_force_comm;
}

void
Space_charge_2d_open_hockney::auto_tune_comm(bool verbose)
{
    bool output = verbose && (comm2.get_rank() == 0);
    double t0, t1;

    std::vector<double > size(3), offset(3);
    size[0] = size[1] = size[2] = 1.0;
    offset[0] = offset[1] = offset[2] = 0.0;

    domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                    size, offset, grid_shape, periodic_z));
    set_fixed_domain(domain_sptr);
    Rectangular_grid fake_local_charge_density(doubled_domain_sptr);
    if (output) {
        std::cout
                << "Space_charge_2d_open_hockney::auto_tune_comm: trying get_global_charge_density2_reduce_scatter\n";
    }
    MPI_Barrier(comm2.get());
    t0 = MPI_Wtime();
    get_global_charge_density2_reduce_scatter(fake_local_charge_density);
    MPI_Barrier(comm2.get());
    t1 = MPI_Wtime();
    if (output) {
        std::cout << "Space_charge_2d_open_hockney::auto_tune_comm: time = "
                << t1 - t0 << " seconds\n";
    }
    double best_time = t1 - t0;
    charge_density_comm = reduce_scatter;

    if (output) {
        std::cout
                << "Space_charge_2d_open_hockney::auto_tune_comm: trying get_global_charge_density2_allreduce\n";
    }
    MPI_Barrier(comm2.get());
    t0 = MPI_Wtime();
    get_global_charge_density2_allreduce(fake_local_charge_density);
    MPI_Barrier(comm2.get());
    t1 = MPI_Wtime();
    if (output) {
        std::cout << "Space_charge_2d_open_hockney::auto_tune_comm: time = "
                << t1 - t0 << " seconds\n";
    }
    if ((t1 - t0) < best_time) {
        charge_density_comm = charge_allreduce;
    }

    Distributed_rectangular_grid fake_local_e_force(doubled_domain_sptr, 
            doubled_lower, doubled_upper, comm2);

#if 1
    if (output) {
        std::cout
                << "Space_charge_2d_open_hockney::auto_tune_comm: trying get_global_electric_force2_gatherv_bcast\n";
    }
    MPI_Barrier(comm2.get());
    t0 = MPI_Wtime();
    get_global_electric_force2_gatherv_bcast(fake_local_e_force);
    MPI_Barrier(comm2.get());
    t1 = MPI_Wtime();
    if (output) {
        std::cout << "Space_charge_2d_open_hockney::auto_tune_comm: time = "
                << t1 - t0 << " seconds\n";
    }
    best_time = t1 - t0;
    e_force_comm = gatherv_bcast;
#endif

    if (output) {
        std::cout
                << "Space_charge_2d_open_hockney::auto_tune_comm: trying get_global_electric_force2_allgatherv\n";
    }
    MPI_Barrier(comm2.get());
    t0 = MPI_Wtime();
    get_global_electric_force2_allgatherv(fake_local_e_force);
    MPI_Barrier(comm2.get());
    t1 = MPI_Wtime();
    if (output) {
        std::cout << "Space_charge_2d_open_hockney::auto_tune_comm: time = "
                << t1 - t0 << " seconds\n";
    }
    if ((t1 - t0) < best_time) {
        best_time = t1 - t0;
        e_force_comm = allgatherv;
    }

    if (output) {
        std::cout
                << "Space_charge_2d_open_hockney::auto_tune_comm: trying get_global_electric_force2_allreduce\n";
    }
    MPI_Barrier(comm2.get());
    t0 = MPI_Wtime();
    get_global_electric_force2_allreduce(fake_local_e_force);
    MPI_Barrier(comm2.get());
    t1 = MPI_Wtime();
    if (output) {
        std::cout << "Space_charge_2d_open_hockney::auto_tune_comm: time = "
                << t1 - t0 << " seconds\n";
    }
    if ((t1 - t0) < best_time) {
        e_force_comm = e_force_allreduce;
    }

    if (output) {
        std::cout
                << "Space_charge_2d_open_hockney::auto_tune_comm: selected charge_density_comm = "
                << charge_density_comm << std::endl
                << "Space_charge_2d_open_hockney::auto_tune_comm: selected e_force_comm = "
                << e_force_comm << std::endl;
    }

    // this would be unset_fixed_domain(), if such a method existed
    domain_fixed = false;
    have_domains = false;
    this->doubled_domain_sptr.reset();
    this->domain_sptr.reset();
}

void
Space_charge_2d_open_hockney::set_doubled_domain()
{
    std::vector<double > doubled_size(3);
    for (int i = 0; i < 2; ++i) {
        doubled_size[i] = 2 * domain_sptr->get_physical_size()[i];
    }
    doubled_size[2] = domain_sptr->get_physical_size()[2];

    doubled_domain_sptr = Rectangular_grid_domain_sptr(
            new Rectangular_grid_domain(doubled_size,
                    domain_sptr->get_physical_offset(), doubled_grid_shape,
                    periodic_z));
}

void
Space_charge_2d_open_hockney::set_fixed_domain(
        Rectangular_grid_domain_sptr domain_sptr)
{
    if ((domain_sptr->get_grid_shape()[0] != grid_shape[0])
            || (domain_sptr->get_grid_shape()[1] != grid_shape[1])
            || (domain_sptr->get_grid_shape()[2] != grid_shape[2])) {
        throw runtime_error(
                "Space_charge_2d_open_hockney::set_fixed_domain requires a shape\nequal to that of the parent object, but with zyx ordering.");
    }
    this->domain_sptr = domain_sptr;
    set_doubled_domain();
    domain_fixed = true;
    have_domains = true;
}

void
Space_charge_2d_open_hockney::update_domain(Bunch const& bunch)
{
    if (!domain_fixed) {
        MArray1d mean(Diagnostics::calculate_mean(bunch));
        MArray1d std(Diagnostics::calculate_std(bunch, mean));
        std::vector<double > size(3);
        std::vector<double > offset(3);
        // domain is in xyz order
        offset[0] = mean[Bunch::x];
        size[0] = n_sigma * std[Bunch::x];
        offset[1] = mean[Bunch::y];
        size[1] = n_sigma * std[Bunch::y];
        if (grid_entire_period) {
            offset[2] = 0.0;
            size[2] = z_period;
        } else {
            offset[2] = mean[Bunch::z];
            size[2] = n_sigma * std[Bunch::z];
        }
        domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                size, offset, grid_shape, periodic_z));
        set_doubled_domain();
        have_domains = true;
    }
}

Rectangular_grid_domain_sptr
Space_charge_2d_open_hockney::get_domain_sptr() const
{
    if (!have_domains) {
        throw runtime_error(
                "Space_charge_2d_open_hockney::get_domain_sptr: domain not set");
    }
    return domain_sptr;
}

Rectangular_grid_domain_sptr
Space_charge_2d_open_hockney::get_doubled_domain_sptr() const
{
    if (!have_domains) {
        throw runtime_error(
                "Space_charge_2d_open_hockney::get_doubled_domain_sptr: domain not set");
    }
    return doubled_domain_sptr;
}

Rectangular_grid_sptr
Space_charge_2d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
    update_domain(bunch);
    bunch_particle_charge = bunch.get_particle_charge();
    bunch_real_num = bunch.get_real_num();
    bunch_total_num = bunch.get_total_num();
    Rectangular_grid_sptr local_rho_sptr(new Rectangular_grid(doubled_domain_sptr));
    deposit_charge_rectangular_2d(*local_rho_sptr, bunch);
    return local_rho_sptr;
}

Distributed_rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_charge_density2_reduce_scatter(
        Rectangular_grid const& local_charge_density)
{
    // jfa: here is where we do something complicated, but (potentially) efficient
    // in calculating a version of the charge density that is just global enough
    // to fill in the doubled global charge density
    // dest_array stores the portion of global charge density needed on each
    // processor. It has to have the same shape in the non-distributed 
    // dimensions as the charge density in order to work with MPI_Reduce_scatter.
    const std::complex<double > * source_2dc
            = local_charge_density.get_grid_points_2dc().origin();
    std::complex<double > * dest_2dc;
    MArray2dc dest_array_2dc(boost::extents[1][1]);
    dest_array_2dc.resize(boost::extents[extent_range(doubled_lower,
            doubled_upper)][doubled_grid_shape[1]]);
    dest_2dc = multi_array_offset(dest_array_2dc, doubled_lower, 0);
    int error_2dc = MPI_Reduce_scatter((void *) source_2dc, (void *) dest_2dc,
            &real_lengths[0], MPI_DOUBLE_COMPLEX, MPI_SUM, comm2.get());

    int error_1d = MPI_Allreduce(MPI_IN_PLACE,
            (void*) local_charge_density.get_grid_points_1d().origin(),
            local_charge_density.get_grid_points_1d().num_elements(),
            MPI_DOUBLE, MPI_SUM, comm2.get());
#if 0
    const double * source_1d 
            = local_charge_density.get_grid_points_1d().origin();
    double * dest_1d;
    MArray1d dest_array_1d(boost::extents[1]);
    dest_array_1d.resize(boost::extents[doubled_grid_shape[2]]);
    dest_1d = multi_array_offset(dest_array_1d, 0);

    int error_1d = MPI_Reduce_scatter((void *) source_1d, (void *) dest_1d,
            &real_lengths_1d[0], MPI_DOUBLE, MPI_SUM, comm2.get());
#endif

    if ((error_2dc != MPI_SUCCESS) || (error_1d != MPI_SUCCESS)) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney::get_global_charge_density2_reduce_scatter");
    }
    Distributed_rectangular_grid_sptr rho2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr,
                    doubled_lower, doubled_upper, doubled_grid_shape, comm2));
    for (int i = rho2->get_lower(); i < rho2->get_upper(); ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            rho2->get_grid_points_2dc()[i][j] = dest_array_2dc[i][j];
        }
    }
    for (int k = 0; k < doubled_grid_shape[2]; ++k) {
        //rho2->get_grid_points_1d()[k] = dest_array_1d[k];
        rho2->get_grid_points_1d()[k] 
                = local_charge_density.get_grid_points_1d()[k];
    }
    rho2->set_normalization(1.0);
    return rho2;
}

Distributed_rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_charge_density2_allreduce(
        Rectangular_grid const& local_charge_density)
{
    int error_2d = MPI_Allreduce(MPI_IN_PLACE,
            (void*) local_charge_density.get_grid_points_2dc().origin(),
            local_charge_density.get_grid_points_2dc().num_elements(), 
            MPI_DOUBLE_COMPLEX, MPI_SUM, comm2.get());
    int error_1d = MPI_Allreduce(MPI_IN_PLACE,
            (void*) local_charge_density.get_grid_points_1d().origin(),
            local_charge_density.get_grid_points_1d().num_elements(),
            MPI_DOUBLE, MPI_SUM, comm2.get());

    if ((error_2d != MPI_SUCCESS) || (error_1d != MPI_SUCCESS)) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney::get_global_charge_density2_allreduce");
    }
    Distributed_rectangular_grid_sptr rho2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr,
                    doubled_lower, doubled_upper, doubled_grid_shape, comm2));
    for (int i = rho2->get_lower(); i < rho2->get_upper(); ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            rho2->get_grid_points_2dc()[i][j]
                    = local_charge_density.get_grid_points_2dc()[i][j];
        }
    }
    for (int k = 0; k < doubled_grid_shape[2]; ++k) {
        rho2->get_grid_points_1d()[k]
                = local_charge_density.get_grid_points_1d()[k];
    }
    rho2->set_normalization(1.0);

    return rho2;
}

Distributed_rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_charge_density2(
        Rectangular_grid const& local_charge_density)
{
    switch (charge_density_comm) {
    case reduce_scatter:
        return get_global_charge_density2_reduce_scatter(local_charge_density);
    case charge_allreduce:
        return get_global_charge_density2_allreduce(local_charge_density);
    default:
        throw runtime_error(
                "Space_charge_2d_open_hockney: invalid charge_density_comm");
    }
}

Distributed_rectangular_grid_sptr
Space_charge_2d_open_hockney::get_green_fn2_pointlike()
{
    if (doubled_domain_sptr == NULL) {
        throw runtime_error(
                "Space_charge_2d_open_hockney::get_green_fn2_pointlike called before domain specified");
    }
    int lower = distributed_fft2d_sptr->get_lower();
    int upper = distributed_fft2d_sptr->get_upper();
    Distributed_rectangular_grid_sptr G2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr, 
                    lower, upper, doubled_grid_shape, comm2));

    double hx = domain_sptr->get_cell_size()[0];
    double hy = domain_sptr->get_cell_size()[1];
    double hz = domain_sptr->get_cell_size()[2];

    const double epsilon = 0.01;
    double dx, dy, Gx, Gy;

    for (int ix = doubled_lower; ix < doubled_upper; ++ix) {
        if (ix > grid_shape[0]) {
            dx = (ix - doubled_grid_shape[0]) * hx;
        } else {
            dx = ix * hx;
        }
        for (int iy = 0; iy < doubled_grid_shape[1]; ++iy) {
            if (iy > grid_shape[1]) {
                dy = (iy - doubled_grid_shape[1]) * hy;
            } else {
                dy = iy * hy;
            }
            if ((ix == grid_shape[0]) || (iy == grid_shape[1])) {
                Gx = 0.0;
                Gy = 0.0;
            } else {
                Gx = dx / (dx * dx + dy * dy + hx * hy * epsilon * epsilon);
                Gy = dy / (dx * dx + dy * dy + hx * hy * epsilon * epsilon);
            }
            G2->get_grid_points_2dc()[ix][iy] = std::complex<double >(Gx, Gy);
        }
    }
    G2->set_normalization(1.0);

    return G2;
}

Distributed_rectangular_grid_sptr
Space_charge_2d_open_hockney::get_local_force2(
        Distributed_rectangular_grid & charge_density2,
        Distributed_rectangular_grid & green_fn2)
{
    double t;
    t = simple_timer_current();

    int lower = distributed_fft2d_sptr->get_lower();
    int upper = distributed_fft2d_sptr->get_upper();

    MArray2dc rho2hat(
            boost::extents[extent_range(lower, upper)][doubled_grid_shape[1]]);
    MArray2dc G2hat(
            boost::extents[extent_range(lower, upper)][doubled_grid_shape[1]]);
    MArray2dc local_force2hat(
            boost::extents[extent_range(lower, upper)][doubled_grid_shape[1]]);

    t = simple_timer_show(t, "sc-fft-setup");
    // FFT
    distributed_fft2d_sptr->transform(charge_density2.get_grid_points_2dc(), 
            rho2hat);
    distributed_fft2d_sptr->transform(green_fn2.get_grid_points_2dc(), G2hat);

    t = simple_timer_show(t, "sc-fft-transform");
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            local_force2hat[i][j] = rho2hat[i][j] * G2hat[i][j];
        }
    }
    Distributed_rectangular_grid_sptr local_force2(
            new Distributed_rectangular_grid(doubled_domain_sptr, 
                    lower, upper, doubled_grid_shape, comm2));
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            local_force2->get_grid_points_2dc()[i][j] = 0.0;
        }
    }
    t = simple_timer_show(t, "sc-fft-multiplication/initialization");
    // inverse FFT
    distributed_fft2d_sptr->inv_transform(local_force2hat, 
            local_force2->get_grid_points_2dc());
    t = simple_timer_show(t, "sc-fft-inv_transform");

    double hx, hy, hz;
    hx = domain_sptr->get_cell_size()[0];
    hy = domain_sptr->get_cell_size()[1];
    hz = domain_sptr->get_cell_size()[2];
    double normalization = hx * hy;  // volume element in integral
    normalization *= 1.0 / (4.0 * pi * epsilon0);
    normalization *= hz;             // dummy factor from weight0 of deposit.cc
    // normalizations for the lamabda factor

    // average line density for real particles
    normalization *= 1.0 / doubled_domain_sptr->get_physical_size()[2];
    // average line density for macro particles
    normalization *= 2.0 * doubled_domain_sptr->get_physical_size()[2]
            / (hz * bunch_total_num);

    normalization *= bunch_particle_charge * pconstants::e;
    normalization *= charge_density2.get_normalization();
    normalization *= green_fn2.get_normalization();
    normalization *= distributed_fft2d_sptr->get_roundtrip_normalization();
    local_force2->set_normalization(normalization);
    t = simple_timer_show(t, "sc-fft-normalization");

    return local_force2;
}

Rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_electric_force2_gatherv_bcast(
        Distributed_rectangular_grid const& dist_force)
{
    Rectangular_grid_sptr global_force2(new Rectangular_grid(
            doubled_domain_sptr));
    const int root = 0;
    int error;
    if (in_group1) {
        int rank = comm2.get_rank();
        error = MPI_Gatherv((void *) (dist_force.get_grid_points_2dc().origin()
                + lowers1[rank]), lengths1[rank], MPI_DOUBLE_COMPLEX,
                (void*) global_force2->get_grid_points_2dc().origin(), 
                &lengths1[0], &lowers1[0], MPI_DOUBLE_COMPLEX, root, 
                comm1.get());
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                    "MPI error in Space_charge_2d_open_hockney(MPI_Gatherv)");
        }
    }
    int total_length = doubled_grid_shape[0] * doubled_grid_shape[1];
    error = MPI_Bcast(global_force2->get_grid_points_2dc().origin(), 
            total_length, MPI_DOUBLE_COMPLEX, root, comm2.get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney(MPI_Bcast)");
    }
    global_force2->set_normalization(dist_force.get_normalization());

    return global_force2;
}

Rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_electric_force2_allgatherv(
        Distributed_rectangular_grid const& dist_force)
{
    Rectangular_grid_sptr global_force2(new Rectangular_grid(
            doubled_domain_sptr));
    std::vector<int > lowers12(comm2.get_size());  // lowers1 on comm2
    std::vector<int > lengths12(comm2.get_size()); // lengths1 on comm2
    int size1 = lowers1.size();
    for (int rank = 0; rank < comm2.get_size(); ++rank) {
        if (rank < size1) {
            lowers12[rank] = lowers1[rank];
            lengths12[rank] = lengths1[rank];
        } else {
            lowers12[rank] = 0;
            lengths12[rank] = 0;
        }
    }
    int rank = comm2.get_rank();
    int error = MPI_Allgatherv((void *) (
            dist_force.get_grid_points_2dc().origin() + lowers12[rank]), 
            lengths12[rank], MPI_DOUBLE_COMPLEX,
            (void*) global_force2->get_grid_points_2dc().origin(), 
            &lengths12[0], &lowers12[0], MPI_DOUBLE_COMPLEX, comm2.get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney(MPI_Allgatherv)");
    }
    global_force2->set_normalization(dist_force.get_normalization());

    return global_force2;
}

Rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_electric_force2_allreduce(
        Distributed_rectangular_grid const& dist_force)
{
    Rectangular_grid_sptr global_force2(new Rectangular_grid(
            doubled_domain_sptr));
    for (int i = 0; i < doubled_grid_shape[0]; ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            global_force2->get_grid_points_2dc()[i][j] = 0.0;
        }
    }
    for (int i = dist_force.get_lower(); i < std::min(doubled_grid_shape[0],
            dist_force.get_upper()); ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            global_force2->get_grid_points_2dc()[i][j]
                    = dist_force.get_grid_points_2dc()[i][j];
        }
    }
    int error = MPI_Allreduce(MPI_IN_PLACE,
            (void*) global_force2->get_grid_points_2dc().origin(),
            global_force2->get_grid_points_2dc().num_elements(), 
            MPI_DOUBLE_COMPLEX, MPI_SUM, comm2.get());
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney(MPI_Allreduce in get_global_electric_force2_allreduce)");
    }
    global_force2->set_normalization(dist_force.get_normalization());

    return global_force2;
}

Rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_electric_force2(
        Distributed_rectangular_grid const& dist_force)
{
    switch (e_force_comm) {
    case gatherv_bcast:
        return get_global_electric_force2_gatherv_bcast(dist_force);
    case allgatherv:
        return get_global_electric_force2_allgatherv(dist_force);
    case e_force_allreduce:
        return get_global_electric_force2_allreduce(dist_force);
    default:
        throw runtime_error(
                "Space_charge_2d_open_hockney: invalid e_force_comm");
    }
}

void
Space_charge_2d_open_hockney::apply_kick(Bunch & bunch, 
        Distributed_rectangular_grid const& rho2, 
        Rectangular_grid const& Fn, double delta_t) 
{
    // $\delta \vec{p} = \vec{F} \delta t = q \vec{E} \delta t$
    double q = bunch.get_particle_charge() * pconstants::e; // [C]
    // delta_t_beam: [s] in beam frame
    double delta_t_beam = delta_t / bunch.get_reference_particle().get_gamma();
    // unit_conversion: [N] = [kg m/s] to [Gev/c]
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
    // scaled p = p/p_ref
    double p_scale = 1.0 / bunch.get_reference_particle().get_momentum();
    double factor = unit_conversion * delta_t_beam * Fn.get_normalization()
            * p_scale;
    Rectangular_grid_domain & domain(*Fn.get_domain_sptr());
    MArray2dc_ref grid_points(Fn.get_grid_points_2dc());
    MArray1d_ref grid_points_1d(rho2.get_grid_points_1d());
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        double x = bunch.get_local_particles()[part][Bunch::x];
        double y = bunch.get_local_particles()[part][Bunch::y];
        double z = bunch.get_local_particles()[part][Bunch::z];
        std::vector<double > grid_val = interpolate_rectangular_2d(x, y, z, 
                domain, grid_points, grid_points_1d); 
        bunch.get_local_particles()[part][1] += factor * grid_val[0];
        bunch.get_local_particles()[part][3] += factor * grid_val[1];
    }
}

void
Space_charge_2d_open_hockney::apply(Bunch & bunch, double time_step,
        Step & step)
{
    double t;
    t = simple_timer_current();
    bunch.convert_to_state(Bunch::fixed_t);
    t = simple_timer_show(t, "sc-convert-to-state");
    Rectangular_grid_sptr local_rho(get_local_charge_density(bunch)); // [C/m^3]
    t = simple_timer_show(t, "sc-get-local-rho");
    Distributed_rectangular_grid_sptr rho2(get_global_charge_density2(
            *local_rho)); // [C/m^3]
    local_rho.reset();
    t = simple_timer_show(t, "sc-get-global-rho");
    Distributed_rectangular_grid_sptr G2(get_green_fn2_pointlike());
    t = simple_timer_show(t, "sc-get-green-fn");
    Distributed_rectangular_grid_sptr local_force2(get_local_force2(*rho2, 
            *G2));        // [N]
    //rho2.reset();
    G2.reset();
    t = simple_timer_show(t, "sc-get-local_force");
    Rectangular_grid_sptr Fn(get_global_electric_force2(*local_force2)); // [N]
    local_force2.reset();
    t = simple_timer_show(t, "sc-get-global-force");
    bunch.periodic_sort(Bunch::z);
    t = simple_timer_show(t, "sc-sort");
    apply_kick(bunch, *rho2, *Fn, time_step);
    rho2.reset();
    Fn.reset();
    t = simple_timer_show(t, "sc-apply-kick");
}

Space_charge_2d_open_hockney::~Space_charge_2d_open_hockney()
{
    int error;
    if (in_group1) {
        error = MPI_Comm_free(&mpi_comm1);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                    "MPI error in Space_charge_2d_open_hockney(MPI_Comm_free)");
        }
    }
    error = MPI_Group_free(&group1);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney(MPI_Group_free(1))");
    }
    error = MPI_Group_free(&group2);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney(MPI_Group_free(2))");
    }

}
