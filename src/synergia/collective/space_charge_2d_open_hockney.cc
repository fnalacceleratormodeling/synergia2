
#include "space_charge_2d_open_hockney.h"

#include "synergia/bunch/core_diagnostics.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/foundation/physical_constants.h"

#include "deposit.h"

//#include "interpolate_rectangular_zyx.h"
//#include "synergia/utils/multi_array_offsets.h"
//#include "synergia/utils/simple_timer.h"
//#include "synergia/utils/hdf5_writer.h"

using mconstants::pi;
using pconstants::epsilon0;

namespace
{
    double
    get_smallest_non_tiny( double val, 
                           double other1, 
                           double other2, 
                           double tiny )
    {
        double retval;
        if (val > tiny) 
        {
            retval = val;
        } 
        else 
        {
            if ((other1 > tiny) && (other2 > tiny)) 
            {
                retval = std::min(other1, other2);
            } 
            else 
            {
                retval = std::max(other1, other2);
            }
        }

        return retval;
    }
}


Space_charge_2d_open_hockney::Space_charge_2d_open_hockney(
        Space_charge_2d_open_hockney_options ops)
    : Collective_operator("sc_2d_open_hockney", 1.0)
    , options(ops)
    , domain(ops.shape, {1.0, 1.0, 1.0})
    , doubled_domain(ops.doubled_shape, {1.0, 1.0, 1.0})
{
}

void 
Space_charge_2d_open_hockney::apply_impl(
            Bunch_simulator & simulator, 
            double time_step, 
            Logger & logger)
{
    logger << "    Space charge 2d open hockney\n";
    apply_bunch(simulator[0][0], time_step, logger);

#if 0
    if (bunch.get_total_num() <= 1) return;

    //setup_communication(bunch.get_comm_sptr());

    double t, t_total;
    t_total = simple_timer_current();
    t = t_total;

    //update_domain(bunch);
    t = simple_timer_show(t, "get_local_rho-update_domain");

    auto local_rho = get_local_charge_density(bunch); // [C/m^3]
    t = simple_timer_show(t, "sc2doh-get_local_rho");

    auto rho2 = get_global_charge_density2(local_rho, bunch.get_comm());
    t = simple_timer_show(t, "sc2doh-get_global_rho");

    auto G2 = get_green_fn2_pointlike();
    t = simple_timer_show(t, "sc2doh-get_green_fn");

    auto local_force2 = get_local_force2(rho2, G2); // [N]
    t = simple_timer_show(t, "sc2doh-get_local_force");

    auto Fn = get_global_electric_force2(local_force2); // [N]
    t = simple_timer_show(t, "sc2doh-get_global_force");

    apply_kick(bunch, rho2, Fn, time_step);
    t = simple_timer_show(t, "sc2doh-apply_kick");

    t_total = simple_timer_show(t_total, "collective_operator_apply-sc2doh");
#endif
}

void
Space_charge_2d_open_hockney::apply_bunch(
            Bunch & bunch, 
            double time_step, 
            Logger & logger)
{

    //if (!bunch.get_comm().has_this_rank()) return;

    //setup_communication(bunch.get_comm());

    update_domain(bunch);

    auto local_rho = get_local_charge_density(bunch); // [C/m^3]

#if 0
    for(int i=24; i<32; ++i)
    {
        logger << local_rho(32*64 + i) << ", ";
    }
    logger << "\n";
#endif
}

void
Space_charge_2d_open_hockney::update_domain(Bunch const & bunch)
{
    auto mean = Core_diagnostics::calculate_mean(bunch);
    auto std  = Core_diagnostics::calculate_std(bunch, mean);

    const double tiny = 1.0e-10;

    if ((std[Bunch::x] < tiny) 
            && (std[Bunch::y] < tiny) 
            && (std[Bunch::z] < tiny)) 
    {
        throw std::runtime_error(
                "Space_charge_3d_open_hockney_eigen::update_domain: "
                "all three spatial dimensions have neglible extent");
    }

    std::array<double, 3> offset { 
        mean[0], 
        mean[2], 
        mean[4] };

    std::array<double, 3> size { 
        options.n_sigma * get_smallest_non_tiny(std[0], std[2], std[4], tiny),
        options.n_sigma * get_smallest_non_tiny(std[2], std[0], std[4], tiny),
        options.n_sigma * get_smallest_non_tiny(std[4], std[0], std[2], tiny) };

    std::array<double, 3> doubled_size { 
        size[0] * 2.0, 
        size[1] * 2.0, 
        size[2] };

    domain = Rectangular_grid_domain(
            options.shape, size, offset, false);

    doubled_domain = Rectangular_grid_domain(
            options.doubled_shape, doubled_size, offset, false);
}


karray1d
Space_charge_2d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
    particle_bin = karray2d_dev("bin", bunch.get_local_num(), 6);
    return deposit_charge_rectangular_2d_kokkos(doubled_domain, particle_bin, bunch);
}

#if 0
Rectangular_grid
Space_charge_2d_open_hockney::get_global_charge_density2(
        Rectangular_grid_2d const & rho_xy, 
        Rectangular_grid_1d const & rho_z, 
        Commxx_sptr comm_sptr)
{
    //setup_communication(comm_sptr);

    int error_2d = MPI_Allreduce(MPI_IN_PLACE,
            (void*) rho_xy.data(), rho_xy.span() * 2,
            MPI_DOUBLE, MPI_SUM, comm_sptr->get());

    int error_1d = MPI_Allreduce(MPI_IN_PLACE,
            (void*)rho_z.data(), rho_z.span(),
            MPI_DOUBLE, MPI_SUM, comm_sptr->get());

    if ((error_2d != MPI_SUCCESS) || (error_1d != MPI_SUCCESS)) 
    {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney::"
                "get_global_charge_density2_allreduce");
    }

    Distributed_rectangular_grid_sptr rho2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr,
                    doubled_lower, doubled_upper, doubled_grid_shape,
                    comm_sptr));

    #pragma omp parallel for
    for (int i = rho2->get_lower(); i < rho2->get_upper(); ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            rho2->get_grid_points_2dc()[i][j]
                    = local_charge_density.get_grid_points_2dc()[i][j];
        }
    }

    #pragma omp parallel for
    for (int k = 0; k < doubled_grid_shape[2]; ++k) {
        rho2->get_grid_points_1d()[k]
                = local_charge_density.get_grid_points_1d()[k];
    }

    rho2->set_normalization(1.0);

    return rho2;
}
#endif


#if 0
void
Space_charge_2d_open_hockney::setup_communication(
        Commxx_sptr const& bunch_comm_sptr)
{
    if (comm2_sptr != commxx_divider_sptr->get_commxx_sptr(bunch_comm_sptr)) {
        comm2_sptr = commxx_divider_sptr->get_commxx_sptr(bunch_comm_sptr);
        setup_derived_communication();
    }
}

void
Space_charge_2d_open_hockney::setup_derived_communication()
{
    distributed_fft2d_sptr = Distributed_fft2d_sptr(
            new Distributed_fft2d(doubled_grid_shape, comm2_sptr));
    std::vector<int > ranks1; // ranks with data from the undoubled domain
    int lower = 0;
    for (int rank = 0; rank < comm2_sptr->get_size(); ++rank) {
        int uppers2 = distributed_fft2d_sptr->get_uppers()[rank];
        int length0;
        if (rank > 0) {
            length0 = uppers2 - distributed_fft2d_sptr->get_uppers()[rank - 1];
        } else {
            length0 = uppers2;
        }
        if (length0 > 0) {
            ranks1.push_back(rank);
            lowers1.push_back(lower);
            int total_length = length0 * doubled_grid_shape[1];
            lengths1.push_back(total_length);
            lower += total_length;
        }
    }
    comm1_sptr = Commxx_sptr(new Commxx(comm2_sptr, ranks1));

    std::vector<int > real_uppers(distributed_fft2d_sptr->get_uppers());
    real_lengths = distributed_fft2d_sptr->get_lengths();
    real_lengths_1d = distributed_fft2d_sptr->get_lengths_1d();
    for (int i = 0; i < comm2_sptr->get_size(); ++i) {
        if (i == 0) {
            real_lengths[0] = real_uppers[0] * doubled_grid_shape[1];
            real_lengths_1d[0] = doubled_grid_shape[2];
        } else {
            real_lengths[i] = (real_uppers[i] - real_uppers[i - 1])
                    * doubled_grid_shape[1];
            real_lengths_1d[i] = doubled_grid_shape[2];
        }
    }
    int my_rank = comm2_sptr->get_rank();
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
}

Space_charge_2d_open_hockney::Space_charge_2d_open_hockney(
        Commxx_divider_sptr commxx_divider_sptr, std::vector<int > const & grid_shape,
        bool need_state_conversion, bool periodic_z, double z_period,
        bool grid_entire_period, double n_sigma) :
                Collective_operator("space charge 2D open hockney", 1.0),
                grid_shape(3),
                doubled_grid_shape(3),
                periodic_z(periodic_z),
                z_period(z_period),
                grid_entire_period(grid_entire_period),
                commxx_divider_sptr(commxx_divider_sptr),
                comm1_sptr(),
                comm2_sptr(),
                n_sigma(n_sigma),
                domain_fixed(false),
                have_domains(false),
                need_state_conversion(need_state_conversion),
                use_cell_coords(true),
                exfile(""),
                eyfile(""),
                dumped(true)
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

    //distributed_fft2d_sptr = Distributed_fft2d_sptr(
    //        new Distributed_fft2d(doubled_grid_shape, comm_sptr));
    //setup_nondoubled_communication();
    setup_default_options();
}

Space_charge_2d_open_hockney::Space_charge_2d_open_hockney(
        Commxx_sptr comm_sptr, std::vector<int > const & grid_shape,
        bool need_state_conversion, bool periodic_z, double z_period,
        bool grid_entire_period, double n_sigma) :
                Collective_operator("space charge 2D open hockney"),
                grid_shape(3),
                doubled_grid_shape(3),
                periodic_z(periodic_z),
                z_period(z_period),
                grid_entire_period(grid_entire_period),
                commxx_divider_sptr(new Commxx_divider),
                comm1_sptr(),
                comm2_sptr(),
                n_sigma(n_sigma),
                domain_fixed(false),
                have_domains(false),
                need_state_conversion(need_state_conversion),
                use_cell_coords(true),
                exfile(""),
                eyfile(""),
                dumped(true)
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
    //distributed_fft2d_sptr = Distributed_fft2d_sptr(
    //        new Distributed_fft2d(doubled_grid_shape, comm_sptr));
    //setup_nondoubled_communication();
    setup_default_options();
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
                "Space_charge_2d_open_hockney::set_fixed_domain requires a shape\nequal to that of the parent object.");
    }
    this->domain_sptr = domain_sptr;
    set_doubled_domain();
    domain_fixed = true;
    have_domains = true;
}

// get_smallest_non_tiny is in space_charge_3d_open_hockney.cc
double
get_smallest_non_tiny(double val, double other1, double other2, double tiny);

void
Space_charge_2d_open_hockney::update_domain(Bunch const& bunch)
{
    setup_communication(bunch.get_comm_sptr());
    if (!domain_fixed) {
        MArray1d mean(Core_diagnostics::calculate_spatial_mean(bunch));
        MArray1d std(Core_diagnostics::calculate_spatial_std(bunch, mean));
        std::vector<double > size(3);
        std::vector<double > offset(3);
        // domain is in xyz order
        const double tiny = 1.0e-10;
        // protect against degenerate domain sizes
        if ((std[0] < tiny) && (std[1] < tiny)
                && (std[2] < tiny)) {
            throw std::runtime_error(
                    "Space_charge_3d_open_hockney::update_domain: all three spatial dimensions have neglible extent");
        }

        for (int i = 0; i < 3; ++i) {
            offset[i] = mean[i];
            size[i] = n_sigma * std[i];
        }
        if (grid_entire_period) {
        	offset[2] = 0.0;
        	size[2] = z_period;
        }  else {
        	offset[2] = mean[2];
        	size[2] = n_sigma
        			* get_smallest_non_tiny(std[2], std[0],
        					std[2], tiny);
        }
        offset[1] = mean[1];
        size[1] = n_sigma
        		* get_smallest_non_tiny(std[1], std[0],
        				std[2], tiny);
        offset[0] = mean[0];
        size[0] = n_sigma
        		* get_smallest_non_tiny(std[0], std[1],
        				std[2], tiny);

        domain_sptr = Rectangular_grid_domain_sptr(new Rectangular_grid_domain(
                size, offset, grid_shape, periodic_z));
        set_doubled_domain();
        have_domains = true;
    }
}


Rectangular_grid_sptr
Space_charge_2d_open_hockney::get_local_charge_density(Bunch const& bunch)
{
    double t = simple_timer_current();

    update_domain(bunch);
    t = simple_timer_show(t, "get_local_rho-update_domain");

    bunch_particle_charge = bunch.get_particle_charge();
    bunch_total_num = bunch.get_total_num();
    beta = bunch.get_reference_particle().get_beta();
    gamma = bunch.get_reference_particle().get_gamma();

    Rectangular_grid_sptr local_rho_sptr(new Rectangular_grid(doubled_domain_sptr));
    t = simple_timer_show(t, "get_local_rho-setup");

    if (!use_cell_coords) 
    {
        deposit_charge_rectangular_2d(*local_rho_sptr, bunch);
    } 
    else 
    {
        particle_bin_sptr
                = boost::shared_ptr<Raw_MArray2d >(
                        new Raw_MArray2d(boost::extents[bunch.get_local_num()][6]));
        deposit_charge_rectangular_2d_omp_reduce(*local_rho_sptr, *particle_bin_sptr,
                bunch);
    }

    t = simple_timer_show(t, "get_local_rho-deposit");

    return local_rho_sptr;
}

Distributed_rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_charge_density2(
        Rectangular_grid const& local_charge_density, 
        Commxx_sptr comm_sptr)
{
    setup_communication(comm_sptr);
    int error_2d = MPI_Allreduce(MPI_IN_PLACE,
            (void*) local_charge_density.get_grid_points_2dc().origin(),
            local_charge_density.get_grid_points_2dc().num_elements() * 2,
            MPI_DOUBLE, MPI_SUM, comm_sptr->get());
    int error_1d = MPI_Allreduce(MPI_IN_PLACE,
            (void*) local_charge_density.get_grid_points_1d().origin(),
            local_charge_density.get_grid_points_1d().num_elements(),
            MPI_DOUBLE, MPI_SUM, comm_sptr->get());

    if ((error_2d != MPI_SUCCESS) || (error_1d != MPI_SUCCESS)) 
    {
        throw std::runtime_error(
                "MPI error in Space_charge_2d_open_hockney::"
                "get_global_charge_density2_allreduce");
    }

    Distributed_rectangular_grid_sptr rho2 = Distributed_rectangular_grid_sptr(
            new Distributed_rectangular_grid(doubled_domain_sptr,
                    doubled_lower, doubled_upper, doubled_grid_shape,
                    comm_sptr));

    #pragma omp parallel for
    for (int i = rho2->get_lower(); i < rho2->get_upper(); ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            rho2->get_grid_points_2dc()[i][j]
                    = local_charge_density.get_grid_points_2dc()[i][j];
        }
    }

    #pragma omp parallel for
    for (int k = 0; k < doubled_grid_shape[2]; ++k) {
        rho2->get_grid_points_1d()[k]
                = local_charge_density.get_grid_points_1d()[k];
    }

    rho2->set_normalization(1.0);

    return rho2;
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
                    lower, upper, doubled_grid_shape, comm2_sptr));

    double hx = domain_sptr->get_cell_size()[0];
    double hy = domain_sptr->get_cell_size()[1];

    const double epsilon = 0.01;
    double dx, dy, Gx, Gy;

    #pragma omp parallel for private( dx, dy, Gx, Gy )
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

    Raw_MArray2dc rho2hat(
            boost::extents[extent_range(lower, upper)][doubled_grid_shape[1]]);
    Raw_MArray2dc G2hat(
            boost::extents[extent_range(lower, upper)][doubled_grid_shape[1]]);
    Raw_MArray2dc local_force2hat(
            boost::extents[extent_range(lower, upper)][doubled_grid_shape[1]]);

    t = simple_timer_show(t, "get_local_force-fft_setup");

    // FFT
    distributed_fft2d_sptr->transform(charge_density2.get_grid_points_2dc(),
            rho2hat.m);
    t = simple_timer_show(t, "get_local_force-get_rhohat");
    distributed_fft2d_sptr->transform(green_fn2.get_grid_points_2dc(), G2hat.m);
    t = simple_timer_show(t, "get_local_force-get_Ghat");

    Distributed_rectangular_grid_sptr local_force2(
            new Distributed_rectangular_grid(doubled_domain_sptr,
                    lower, upper, doubled_grid_shape, comm2_sptr));
    t = simple_timer_show(t, "get_local_force-construct_local_force2");

    #pragma omp parallel for
    for (int i = lower; i < upper; ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            local_force2hat.m[i][j] = rho2hat.m[i][j] * G2hat.m[i][j];
        }
    }
    t = simple_timer_show(t, "get_local_force-convolution");

    // inverse FFT
    distributed_fft2d_sptr->inv_transform(local_force2hat.m,
            local_force2->get_grid_points_2dc());
    t = simple_timer_show(t, "get_local_force-get_local_force2");

    double hx, hy, hz;
    hx = domain_sptr->get_cell_size()[0];
    hy = domain_sptr->get_cell_size()[1];
    hz = domain_sptr->get_cell_size()[2];
    double normalization = hx * hy;  // volume element in integral
    normalization *= 1.0 / (4.0 * pi * epsilon0);
    normalization *= hz;             // dummy factor from weight0 of deposit.cc
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
    t = simple_timer_show(t, "get_local_force-get_force");

    return local_force2;
}


Rectangular_grid_sptr
Space_charge_2d_open_hockney::get_global_electric_force2_allreduce(
        Distributed_rectangular_grid const& dist_force)
{
    Rectangular_grid_sptr global_force2(new Rectangular_grid(
            doubled_domain_sptr));

    std::memset( (void*)global_force2->get_grid_points_2dc().data(), 0,
            global_force2->get_grid_points_2dc().num_elements()*sizeof(double) );

    #pragma omp parallel for
    for (int i = dist_force.get_lower(); i < std::min(doubled_grid_shape[0],
            dist_force.get_upper()); ++i) {
        for (int j = 0; j < doubled_grid_shape[1]; ++j) {
            global_force2->get_grid_points_2dc()[i][j]
                    = dist_force.get_grid_points_2dc()[i][j];
        }
    }

    int error = MPI_Allreduce(MPI_IN_PLACE,
            (void*) global_force2->get_grid_points_2dc().origin(),
            global_force2->get_grid_points_2dc().num_elements() * 2,
            MPI_DOUBLE, MPI_SUM, comm2_sptr->get());
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
    Rectangular_grid_sptr global_electric_field_sptr;

    switch (e_force_comm) {
    case gatherv_bcast:
        global_electric_field_sptr = get_global_electric_force2_gatherv_bcast(dist_force);
        break;
    case allgatherv:
        global_electric_field_sptr = get_global_electric_force2_allgatherv(dist_force);
        break;
    case e_force_allreduce:
        global_electric_field_sptr = get_global_electric_force2_allreduce(dist_force);
        break;
    default:
        throw runtime_error(
                "Space_charge_2d_open_hockney: invalid e_force_comm");
    }

    // make sure the field is zero at the edge of the grid
    MArray2dc_ref grid_points = global_electric_field_sptr->get_grid_points_2dc();
    int g0 = grid_points.shape()[0];
    int g1 = grid_points.shape()[1];

    for (int i=0; i<g0; ++i)
    {
        grid_points[i][0] = 0.0;
        grid_points[i][g1-1] = 0.0;
    }

    for (int j=0; j<g1; ++j)
    {
        grid_points[0][j] = 0.0;
        grid_points[g0-1][j] = 0.0;
    }

    return global_electric_field_sptr;
}

void
Space_charge_2d_open_hockney::apply_kick(Bunch & bunch,
        Distributed_rectangular_grid const& rho2,
        Rectangular_grid const& Fn, double delta_t)
{
    // EGS ported AM changes for kicks in lab frame from 3D solver
    //AM: kicks are done in the z_lab frame
    // $\delta \vec{p} = \vec{F} \delta t = q \vec{E} \delta t$
    // delta_t_beam: [s] in beam frame
    //  See chapter 11, jackson electrodynamics, for field transformation from bunch frame (BF)
    //  to the lab frame (LF). Keep in mind that \vec{B}_BF=0.
    //  Ex_LF=gamma*Ex_BF, Ey_LF=gamma*Ey_BF, Ez_LF=Ez_BF
    //  Bx_LF=gamma*beta*Ey_BF, By_LF=-gamma*beta*Ex_BF, Bz_LF=Bz_BF=0
    //  Transverse Lorentz force in the lab frame: Fx_LF=q*(Ex_L-beta_z*By_LF)=q*gamma*(1-beta*beta_z)*Ex_BF
    //  Longitudinal Lorentz force in the lab frame:
    //        Fz=q*(Ez_LF+beta_x*By_LF-beta_y*Bx_LF)=q*(Ez_BF-gamma*beta*(beta_x*Ex_BF+beta_y*Ey_BF ))
    // In order to get a conservative approximation!:
    // The following approximations are done: beta_z=beta, beta_x=beta_y=0, thus suppresing
    // the particles' movement relative to the reference particle. The same approximation was employed when
    // the field in the bunch frame was calculated.
    // Thus: Fx_LF=q*Ex_BF/gamma, Fz=q*Ez_BF
    double delta_t_beam = delta_t / bunch.get_reference_particle().get_gamma();
    // unit_conversion: [N] = [kg m/s^2] to [Gev/c]
    double unit_conversion = pconstants::c / (1.0e9 * pconstants::e);
    // scaled p = p/p_ref
    double gamma=bunch.get_reference_particle().get_gamma();
    double beta=bunch.get_reference_particle().get_beta();
    double p_scale = 1.0 / bunch.get_reference_particle().get_momentum();
    // gamma*beta factor introduced here when we are no longer going to the t_bunch frame.
    // That factor was introduced in the fixed_z_lab to fixed_t_bunch conversion.  Physically,
    // gamma factor comes from the lorentz expansion longitudinally in the bunch
    // frame and beta comes because the stored coordinate is c*dt whereas the actual
    // domain is beta*c*dt.
    double factor = unit_conversion * delta_t_beam * Fn.get_normalization()
            * p_scale / (gamma * beta);
    Rectangular_grid_domain & domain(*Fn.get_domain_sptr());
    MArray2dc_ref grid_points(Fn.get_grid_points_2dc());
    MArray1d_ref grid_points_1d(rho2.get_grid_points_1d());
    double bin[6];

    #pragma omp parallel for private(bin)
    for (int part = 0; part < bunch.get_local_num(); ++part) {
        std::complex<double > grid_val;
        if (!use_cell_coords) {
            double x = bunch.get_local_particles()[part][Bunch::x];
            double y = bunch.get_local_particles()[part][Bunch::y];
            double z = bunch.get_local_particles()[part][Bunch::z];
            grid_val = interpolate_rectangular_2d(x, y, z, domain, grid_points,
                            grid_points_1d);
        } else {
            for (int i = 0; i < 6; ++i) {
                bin[i] = (*particle_bin_sptr).m[part][i];
            }
            grid_val = interpolate_rectangular_2d(bin, periodic_z, grid_points,
                            grid_points_1d);
        }
        bunch.get_local_particles()[part][1] += factor * grid_val.real();
        bunch.get_local_particles()[part][3] += factor * grid_val.imag();
    }

    // release the particle_bin buffer
    particle_bin_sptr.reset();

    // kick the spectator particles
    #pragma omp parallel for
    for (int part = 0; part < bunch.get_local_spectator_num(); ++part) 
    {
        double x = bunch.get_local_spectator_particles()[part][Bunch::x];
        double y = bunch.get_local_spectator_particles()[part][Bunch::y];
        double z = bunch.get_local_spectator_particles()[part][Bunch::z];

        std::complex<double> grid_val = 
            interpolate_rectangular_2d(x, y, z, domain, grid_points, grid_points_1d);

        bunch.get_local_spectator_particles()[part][1] += factor * grid_val.real();
        bunch.get_local_spectator_particles()[part][3] += factor * grid_val.imag();
    }
}

void
Space_charge_2d_open_hockney::apply_impl(
        Bunch_simulator & simulator,
        double time_step,
        Logger & logger )
{
    if (bunch.get_total_num() > 1) 
    {
        setup_communication(bunch.get_comm_sptr());

        int comm_compare;
        MPI_Comm_compare(comm2_sptr->get(), bunch.get_comm().get(), &comm_compare);

        if ((comm_compare == MPI_UNEQUAL)
                && (charge_density_comm != charge_allreduce)) {
            throw std::runtime_error(
                    "Space_charge_2d_open_hockney: set_charge_density_comm(charge_allreduce) required when comm != bunch comm");
        }

        double t, t_total;
        t_total = simple_timer_current();
        t = t_total;

        Rectangular_grid_sptr local_rho(get_local_charge_density(bunch)); // [C/m^3]
        t = simple_timer_show(t, "sc2doh-get_local_rho");

        Distributed_rectangular_grid_sptr rho2(
                get_global_charge_density2(*local_rho, bunch.get_comm_sptr())); // [C/m^3]

        local_rho.reset();
        t = simple_timer_show(t, "sc2doh-get_global_rho");

        Distributed_rectangular_grid_sptr G2(get_green_fn2_pointlike());
        t = simple_timer_show(t, "sc2doh-get_green_fn");

        Distributed_rectangular_grid_sptr local_force2(get_local_force2(*rho2, *G2));        // [N]
        G2.reset();
        t = simple_timer_show(t, "sc2doh-get_local_force");

        Rectangular_grid_sptr Fn(get_global_electric_force2(*local_force2)); // [N]
        local_force2.reset();
        t = simple_timer_show(t, "sc2doh-get_global_force");

        apply_kick(bunch, *rho2, *Fn, time_step);
        rho2.reset();
        t = simple_timer_show(t, "sc2doh-apply_kick");

        Fn.reset();
        t = simple_timer_show(t, "sc2doh-finalize");

        t_total = simple_timer_show(t_total, "collective_operator_apply-sc2doh");
    }
}
#endif


