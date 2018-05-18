#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/floating_point.h"
#include "synergia/utils/serialization.h"

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/bunch_simulator.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/propagator.h"

#include "synergia/bunch/bunch.h"

#include "synergia/utils/boost_test_mpi_fixture.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture); // needed to initialize MPI

const double tolerance = 1.0e-12;

static std::string elens_channel =
  "beam, particle=proton, energy=0.25+pmass;\
   lens: elens, l=2.0, current=1.2, eenergy=0.01, radius=0.001, gaussian;\
   channel: sequence, l=2.0, refer=entry;\
   lens, at=0.0;\
   endsequence;";

BOOST_AUTO_TEST_CASE(empty_test)
{
    BOOST_CHECK(true);
}


#if 0
BOOST_AUTO_TEST_CASE(elens_propagation)
{
    // read and lattice from string
    MadX_reader reader;
    reader.parse(elens_channel);
    Lattice_sptr lattice_sptr(reader.get_lattice_sptr("channel"));

    // extract the parameters of the electron lens for later calculation
    Lattice_elements lelems(lattice_sptr->get_elements());
    double current;
    double eenergy;
    double radius;
    double len;

    for (Lattice_elements::const_iterator lit=lelems.begin(); lit!= lelems.end(); ++lit) 
    {
        if ((*lit)->get_type() == "elens") 
        {
            (*lit)->set_string_attribute("extractor_type", "libff");

            current = (*lit)->get_double_attribute("current");
            radius = (*lit)->get_double_attribute("radius");
            eenergy = (*lit)->get_double_attribute("eenergy")*1.0e-3;
            len = (*lit)->get_double_attribute("l");
        }
    }

    if (current == 0.0) 
    {
        throw std::runtime_error("elens element parameters not set");
    }

    Reference_particle refpart(lattice_sptr->get_reference_particle());

    const int total_num = 32;
    const double real_num = 0.5e10;
    Commxx_sptr comm_sptr(new Commxx());

    Bunch_sptr bunch_sptr(new Bunch(refpart, total_num, real_num, comm_sptr));

    const double s3o2 = std::sqrt(3.0)/2.0;
    const double s2o2 = std::sqrt(2.0)/2.0;
    const double xoffset = 0.001;
    const double smallfactor = 0.9e-7;
    const double smalloffset = xoffset*smallfactor;

    MArray2d_ref local_particles(bunch_sptr->get_local_particles());

    local_particles[0][0] = xoffset;
    local_particles[1][0] = xoffset*s3o2;
    local_particles[1][2] = xoffset*0.5;
    local_particles[2][0] = xoffset*s2o2;
    local_particles[2][2] = xoffset*s2o2;
    local_particles[3][0] = xoffset*0.5;
    local_particles[3][2] = xoffset*s3o2;

    for (int i=0; i<4; ++i) 
    {
        local_particles[i+4][0] = -local_particles[i][2];
        local_particles[i+4][2] =  local_particles[i][0];
        local_particles[i+8][0] = -local_particles[i][0];
        local_particles[i+8][2] = -local_particles[i][2];
        local_particles[i+12][0] = -local_particles[i+4][0];
        local_particles[i+12][2] = -local_particles[i+4][2];
    }

    // second set of 16 particles is just the like the first 16 but scaled
    // to exercise the small radius calculation
    for (int i=0; i<16; ++i) 
    {
        for (int j=0; j<4; ++j) 
        {
            local_particles[i+16][j] = local_particles[i][j]*smallfactor;
        }
    }

    // set dp/p of all particles to 0.001
    for (int i=0; i<total_num; ++i) 
    {
        local_particles[i][5] = 0.001;
    }

    std::cout << "0\n";
    Bunch_simulator bunch_simulator(bunch_sptr);
    Independent_stepper_sptr stepper(new Independent_stepper(lattice_sptr, 1, 1));

    std::cout << "1\n";
    Propagator propagator(stepper);
    propagator.set_final_checkpoint(false);

    std::cout << "2\n";
    propagator.propagate(bunch_simulator, 1, 1, 1);

    std::cout << "3\n";
    // expected total kick
    double beta_b = refpart.get_beta();
    double gamma_b = refpart.get_gamma();
    double gamma_e = (eenergy + pconstants::me)/pconstants::me;
    double beta_e = std::sqrt(1.0 - 1.0/(gamma_e*gamma_e));
    double betagamma_p = 1.001*refpart.get_beta()*refpart.get_gamma();
    double gamma_p = std::sqrt(betagamma_p*betagamma_p + 1.0);
    double beta_p = betagamma_p/gamma_p;

    double kick = -2.0 * current * len * pconstants::rp *
            (1.0 + beta_e*beta_p) *
            (1.0 - std::exp(-0.5*xoffset*xoffset/(radius*radius)))/
            (pconstants::e*beta_e*beta_p*beta_b*gamma_b*pconstants::c*xoffset);

    double smallkick = -2.0 * current * len * pconstants::rp *
            (1.0 + beta_e*beta_p) *
            (1.0 - std::exp(-0.5*smalloffset*smalloffset/(radius*radius)))/
            (pconstants::e*beta_e*beta_p*beta_b*gamma_b*pconstants::c*smalloffset);

    const double tolerance= 1.0e-10;

    for (int i=0; i<total_num; ++i) 
    {
        double angle = std::atan2(local_particles[i][2], local_particles[i][0]);
        double kick_angle = std::atan2(local_particles[i][3], local_particles[i][1]);

        // check kick angle is opposite of position
        double anglediff = std::abs((angle-kick_angle))/mconstants::pi;
        BOOST_CHECK(floating_point_equal(anglediff, 1.0, tolerance));

        // check size of kick
        double xp = local_particles[i][1];
        double yp = local_particles[i][3];
        double tot_kick = std::sqrt(xp*xp + yp*yp);

        if (i < 16) 
        {
            BOOST_CHECK(floating_point_equal(kick, -tot_kick, tolerance));
        } 
        else 
        {
            BOOST_CHECK(floating_point_equal(smallkick, -tot_kick, tolerance));
        }
    }


}
#endif


