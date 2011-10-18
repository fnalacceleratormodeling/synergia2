#ifndef IMPEDANCE_H_
#define  IMPEDANCE_H_
#include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"

class Impedance : public Collective_operator
{
private:
    std::string pipe_symmetry;
    int z_grid;
    int nstored_turns;
    double orbit_length;
    double wake_factor;
    double bunch_spacing;
    std::string wake_file;
    std::vector<double> z_coord;
    std::vector<double>  x_wake;
    std::vector<double> y_wake;
    std::vector<double> z_wake;

public:

    Impedance(std::string const & wake_file, double const & orbit_length, double const & bunchsp, int const  & zgrid, std::string const & pipe_symmetry, int const nstored_turns);


    int get_z_grid() const;
    double get_orbit_length() const;
    double get_wake_factor() const;
    double get_bunch_spacing() const;
    std::string get_pipe_symmetry() const;
    std::string get_wake_file_name() const;
    std::vector<double> get_z_coord() const;
    std::vector<double> get_x_wake() const;
    std::vector<double> get_y_wake() const;
    std::vector<double> get_z_wake() const;
    virtual int get_nstored_turns() const;

   /*
    get_n_sigma() const;
    void
    set_fixed_domain(Rectangular_grid_domain_sptr domain_sptr);
    void
    update_domain(Bunch const& bunch);
    Rectangular_grid_domain_sptr
    get_domain_sptr() const;
    Rectangular_grid_domain_sptr
    get_doubled_domain_sptr() const;
    /// Returns local charge density on original grid in [C/m^3]
    Rectangular_grid_sptr
    get_local_charge_density(Bunch const& bunch);
    /// Returns global charge density on doubled grid in [C/m^3]
    Distributed_rectangular_grid_sptr
    get_global_charge_density2(Rectangular_grid const& local_charge_density);
    /// Returns Green function on the doubled grid in [1/m^3]
    Distributed_rectangular_grid_sptr
    get_green_fn2_pointlike();
    Distributed_rectangular_grid_sptr
    get_green_fn2_linear();
    Distributed_rectangular_grid_sptr
    get_scalar_field2(Distributed_rectangular_grid & charge_density22,
            Distributed_rectangular_grid & green_fn2);
    Distributed_rectangular_grid_sptr
    extract_scalar_field(Distributed_rectangular_grid const& scalar_field2);
    Distributed_rectangular_grid_sptr
    get_electric_field_component(
            Distributed_rectangular_grid const& scalar_field, int component);
    Rectangular_grid_sptr
    get_global_electric_field_component(
            Distributed_rectangular_grid const& dist_field); */
   // void impedance_kick(Bunch & bunch, double delta_tau,);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step);
    virtual
    ~Impedance();
};

typedef boost::shared_ptr<Impedance> Impedance_sptr;

#endif /* IMPEDANCE_H_ */

