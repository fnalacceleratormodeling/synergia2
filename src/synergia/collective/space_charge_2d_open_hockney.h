#ifndef SPACE_CHARGE_2D_OPEN_HOCKNEY_H_
#define SPACE_CHARGE_2D_OPEN_HOCKNEY_H_
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/distributed_rectangular_grid.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/distributed_fft2d.h"

/// Note: internal grid is stored in [z][y][x] order, but
/// grid shape expects [x][y][z] order.
class Space_charge_2d_open_hockney : public Collective_operator
{
public:
    enum Green_fn_type
    {
        pointlike = 1, bruteforce = 2
    };
    enum Charge_density_comm
    {
        reduce_scatter = 1, charge_allreduce = 2
    };
    enum E_force_comm
    {
        gatherv_bcast = 1, allgatherv = 2, e_force_allreduce = 3
    };
private:
    std::vector<int > grid_shape, doubled_grid_shape;
    Rectangular_grid_domain_sptr domain_sptr, doubled_domain_sptr;
    bool periodic_z;
    double z_period;
    bool grid_entire_period;
    Green_fn_type green_fn_type;
    Charge_density_comm charge_density_comm;
    E_force_comm e_force_comm;
    Distributed_fft2d_sptr distributed_fft2d_sptr;
    Commxx_sptr comm2_sptr, comm1_sptr;
    std::vector<int > lowers1, lengths1;
    int real_lower, real_upper, real_length;
    std::vector<int > real_lengths;
    std::vector<int > real_lengths_1d;
    int doubled_lower, doubled_upper;
    int real_doubled_lower, real_doubled_upper;
    double n_sigma;
    double bunch_particle_charge, bunch_total_num;
    double beta, gamma;
    bool use_cell_coords;
    bool need_state_conversion;
    bool domain_fixed;
    bool have_domains;
    std::string exfile, eyfile;
    bool dumped;
    void
    setup_nondoubled_communication();
    void
    setup_default_options();
    void
    set_doubled_domain();
    boost::shared_ptr<Raw_MArray2d > particle_bin_sptr;

public:
    Space_charge_2d_open_hockney(Commxx_sptr comm_sptr,
            std::vector<int > const & grid_shape,
            bool need_state_conversion = true, bool periodic_z = false,
            double z_period = 0.0, bool grid_entire_period = false,
            double n_sigma = 8.0);
    /// Note: Use Space_charge_2d_open_hockney::get_internal_grid_shape for
    /// Distributed_fft2d.
    Space_charge_2d_open_hockney(Distributed_fft2d_sptr distributed_fft2d_sptr,
            bool need_state_conversion = true, bool periodic_z = false,
            double z_period = 0.0, bool grid_entire_period = false,
            double n_sigma = 8.0);
    Space_charge_2d_open_hockney();
    virtual Space_charge_2d_open_hockney *
    clone();
    bool
    get_need_state_conversion();
    double
    get_n_sigma() const;
    void
    set_green_fn_type(Green_fn_type green_fn_type);
    Green_fn_type
    get_green_fn_type() const;
    void
    set_charge_density_comm(Charge_density_comm charge_density_comm);
    Charge_density_comm
    get_charge_density_comm() const;
    void
    set_e_force_comm(E_force_comm e_force_comm);
    E_force_comm
    get_e_force_comm() const;
    void
    auto_tune_comm(bool verbose = false);
    void
    set_fixed_domain(Rectangular_grid_domain_sptr domain_sptr);
    void
    update_domain(Bunch const& bunch);
    Rectangular_grid_domain const&
    get_domain() const
    {
        return *domain_sptr;
    }
    Rectangular_grid_domain_sptr
    get_domain_sptr() const;
    Rectangular_grid_domain_sptr
    get_doubled_domain_sptr() const;
    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2_reduce_scatter(
            Rectangular_grid const& local_charge_density, Commxx_sptr comm_sptr);
    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2_allreduce(
            Rectangular_grid const& local_charge_density, Commxx_sptr comm_sptr);
    /// Returns local charge density on doubled grid in [C/m^3]
    Rectangular_grid_sptr
    get_local_charge_density(Bunch const& bunch);
    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2(Rectangular_grid const& local_charge_density,
            Commxx_sptr comm_sptr);
    /// Returns Green function on the doubled grid in [1/m^3]
    Distributed_rectangular_grid_sptr
    get_green_fn2_pointlike();
    Distributed_rectangular_grid_sptr
    get_green_fn2_brute_force();
    Distributed_rectangular_grid_sptr
    get_local_force2(Distributed_rectangular_grid & charge_density2,
            Distributed_rectangular_grid & green_fn2);
    Rectangular_grid_sptr
    get_global_electric_force2_gatherv_bcast(
            Distributed_rectangular_grid const& dist_force);
    Rectangular_grid_sptr
    get_global_electric_force2_allgatherv(
            Distributed_rectangular_grid const& dist_force);
    Rectangular_grid_sptr
    get_global_electric_force2_allreduce(
            Distributed_rectangular_grid const& dist_force);
    Rectangular_grid_sptr
    get_global_electric_force2(
            Distributed_rectangular_grid const& dist_force);
    void
    apply_kick(Bunch & bunch, Distributed_rectangular_grid const& rho2_1d,
            Rectangular_grid const& Fn, double delta_tau);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
    void
    set_files(std::string const& xfile, std::string const& yfile);
    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
            ar & BOOST_SERIALIZATION_NVP(comm2_sptr)
                    & BOOST_SERIALIZATION_NVP(grid_shape)
                    & BOOST_SERIALIZATION_NVP(doubled_grid_shape)
                    & BOOST_SERIALIZATION_NVP(use_cell_coords)
                    & BOOST_SERIALIZATION_NVP(need_state_conversion)
                    & BOOST_SERIALIZATION_NVP(periodic_z)
                    & BOOST_SERIALIZATION_NVP(grid_entire_period)
                    & BOOST_SERIALIZATION_NVP(n_sigma)
                    & BOOST_SERIALIZATION_NVP(domain_fixed)
                    & BOOST_SERIALIZATION_NVP(have_domains)
                    & BOOST_SERIALIZATION_NVP(green_fn_type)
                    & BOOST_SERIALIZATION_NVP(charge_density_comm)
                    & BOOST_SERIALIZATION_NVP(e_force_comm)
                    & BOOST_SERIALIZATION_NVP(exfile)
                    & BOOST_SERIALIZATION_NVP(eyfile)
                    & BOOST_SERIALIZATION_NVP(dumped);
        }
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
            ar & BOOST_SERIALIZATION_NVP(comm2_sptr)
                    & BOOST_SERIALIZATION_NVP(grid_shape)
                    & BOOST_SERIALIZATION_NVP(doubled_grid_shape)
                    & BOOST_SERIALIZATION_NVP(use_cell_coords)
                    & BOOST_SERIALIZATION_NVP(need_state_conversion)
                    & BOOST_SERIALIZATION_NVP(periodic_z)
                    & BOOST_SERIALIZATION_NVP(grid_entire_period)
                    & BOOST_SERIALIZATION_NVP(n_sigma)
                    & BOOST_SERIALIZATION_NVP(domain_fixed)
                    & BOOST_SERIALIZATION_NVP(have_domains)
                    & BOOST_SERIALIZATION_NVP(green_fn_type)
                    & BOOST_SERIALIZATION_NVP(charge_density_comm)
                    & BOOST_SERIALIZATION_NVP(e_force_comm)
                    & BOOST_SERIALIZATION_NVP(exfile)
                    & BOOST_SERIALIZATION_NVP(eyfile)
                    & BOOST_SERIALIZATION_NVP(dumped);
            distributed_fft2d_sptr = Distributed_fft2d_sptr(
                    new Distributed_fft2d(doubled_grid_shape, comm2_sptr));
            setup_nondoubled_communication();
        }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    virtual
    ~Space_charge_2d_open_hockney();
};
BOOST_CLASS_EXPORT_KEY(Space_charge_2d_open_hockney)

typedef boost::shared_ptr<Space_charge_2d_open_hockney >
        Space_charge_2d_open_hockney_sptr; // syndoc:include

#endif /* SPACE_CHARGE_2D_OPEN_HOCKNEY_H_ */
