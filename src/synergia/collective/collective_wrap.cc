#include "space_charge_3d_open_hockney.h"
#include <boost/python.hpp>
#include "synergia/utils/container_conversions.h"

using namespace boost::python;

//class Space_charge_3d_open_hockney : public Collective_operator
//{
//
//public:
//    Space_charge_3d_open_hockney(std::vector<int > const & grid_shape,
//            bool periodic_z, Commxx const& comm, double z_period = 0.0,
//            double n_sigma = 8.0);
//    /// Note: Use Space_charge_3d_open_hockney::get_internal_grid_shape for
//    /// Distributed_fft3d.
//    Space_charge_3d_open_hockney(bool periodic_z,
//            Distributed_fft3d_sptr distributed_fft3d_sptr, double z_period =
//                    0.0, double n_sigma = 8.0);
//    double
//    get_n_sigma() const;
//    void
//    set_fixed_domain(Rectangular_grid_domain_sptr domain_sptr);
//    void
//    update_domain(Bunch const& bunch);
//    Rectangular_grid_domain_sptr
//    get_domain_sptr() const;
//    Rectangular_grid_domain_sptr
//    get_doubled_domain_sptr() const;
//    /// Returns local charge density on original grid in [C/m^3]
//    Rectangular_grid_sptr
//    get_local_charge_density(Bunch const& bunch);
//    /// Returns global charge density on doubled grid in [C/m^3]
//    Distributed_rectangular_grid_sptr
//    get_global_charge_density2(Rectangular_grid const& local_charge_density);
//    /// Returns Green function on the doubled grid in [1/m^3]
//    Distributed_rectangular_grid_sptr
//    get_green_fn2();
//    Distributed_rectangular_grid_sptr
//    get_scalar_field2(Distributed_rectangular_grid & charge_density22,
//            Distributed_rectangular_grid & green_fn2);
//    Distributed_rectangular_grid_sptr
//    extract_scalar_field(Distributed_rectangular_grid const& scalar_field2);
//    Distributed_rectangular_grid_sptr
//    get_electric_field_component(
//            Distributed_rectangular_grid const& scalar_field, int component);
//    Rectangular_grid_sptr
//    get_global_electric_field_component(
//            Distributed_rectangular_grid const& dist_field);
//    void
//    apply_kick(Bunch & bunch, Rectangular_grid const& En, double delta_tau,
//            int component);
//    virtual void
//    apply(Bunch & bunch, double time_step, Step & step);
//    ~Space_charge_3d_open_hockney();
//};

BOOST_PYTHON_MODULE(collective)
{
    class_<Space_charge_3d_open_hockney, Space_charge_3d_open_hockney_sptr,
        bases<Collective_operator > >("Space_charge_3d_open_hockney",
                init<std::vector<int > const &, bool, Commxx const&>())
//    class_<Space_charge_3d_open_hockney, Space_charge_3d_open_hockney_sptr>
//        ("Space_charge_3d_open_hockney",
//            init<std::vector<int > const &, bool, Commxx const&>())
        .def("apply", &Space_charge_3d_open_hockney::apply)
        ;
}
