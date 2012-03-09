#ifndef SPACE_CHARGE_RECTANGULAR_H_
#define SPACE_CHARGE_RECTANGULAR_H_
// #include "synergia/simulation/operator.h"
#include "synergia/simulation/step.h"
#include "synergia/bunch/bunch.h"
#include "synergia/collective/rectangular_grid_domain.h"
#include "synergia/collective/rectangular_grid.h"
#include "synergia/collective/distributed_rectangular_grid.h"
#include "synergia/utils/commxx.h"
#include  "synergia/collective/fftw_rectangular_helper.h"

// #include "synergia/utils/distributed_fft3d.h"


class Space_charge_rectangular : public Collective_operator
{
private:
    std::vector<double > pipe_size; //pipe size, x,y,x meters
    std::vector<int > grid_shape;
    Rectangular_grid_domain_sptr domain_sptr;
    Commxx comm_f;
    Fftw_rectangular_helper_sptr fftw_helper_sptr;
    bool have_fftw_helper;
    void
    fill_guards_pplanes(Distributed_rectangular_grid & phi, int lower, int upper, int lengthx,
                          MArray2d & g_lower, MArray2d &g_upper);
    void
    construct_fftw_helper(Commxx const & comm);

public:

    Space_charge_rectangular(Commxx const & comm_f, std::vector<double > const & pipe_size, std::vector<int > const & grid_shape);
    Space_charge_rectangular(std::vector<double > const & pipe_size, std::vector<int > const & grid_shape);

    void
    set_fftw_helper(Commxx const & comm);

    std::vector<double >
    get_pipe_size() const;

    std::vector<int >
    get_grid_shape() const;

   Rectangular_grid_domain_sptr
   get_domain_sptr() const;

   Fftw_rectangular_helper_sptr
   get_fftw_helper_sptr() const;


   Rectangular_grid_sptr
   get_charge_density(Bunch const& bunch);


   Distributed_rectangular_grid_sptr
   get_phi_local(Rectangular_grid & rho, Bunch const& bunch);

   Rectangular_grid_sptr
   get_En( Distributed_rectangular_grid & phi_local,Bunch const& bunch, int component);

    void
    apply_kick(Bunch & bunch, Rectangular_grid const& En, double time_step, int component);

    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
    virtual
    ~Space_charge_rectangular();
};

typedef boost::shared_ptr<Space_charge_rectangular> Space_charge_rectangular_sptr;

#endif /* SPACE_CHARGE_RECTANGULAR_H_ */
