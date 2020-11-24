
#include "synergia/foundation/trigon.h"
#include "synergia/bunch/bunch.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/libFF/ff_element.h"
#include "synergia/utils/lsexpr.h"
#include "synergia/utils/simple_timer.h"


void run()
{
    Logger screen(0, LoggerV::DEBUG);

    auto lsexpr = read_lsexpr_file("sis18-6.lsx");
    Lattice lattice(lsexpr);

    // get the reference particle
    auto const& ref = lattice.get_reference_particle();

    screen << "reference momentum = " 
        << ref.get_momentum() << " GeV\n";

    // closed orbit
    auto probe = Lattice_simulator::calculate_closed_orbit(lattice);

    // trigon bunch
    using trigon_t = Trigon<double, 1, 6>;

    bunch_t<trigon_t> tb(ref);
    auto parts = tb.get_host_particles();

    // init value
    for(int i=0; i<6; ++i) 
        parts(0, i) = trigon_t(probe[i], i);

    // check in
    tb.checkin_particles();

    // propagate trigon
    simple_timer_start("trigon propagate");

    for(auto & ele : lattice.get_elements())
    {
        if (ele.get_type() == element_type::rfcavity)
        {
            ele.set_double_attribute("volt", 0.0);
        }
        else if (ele.get_type() == element_type::quadrupole)
        {
            //ele.set_double_attribute("k1", 0.0);
            //ele.set_double_attribute("k1s", 0.0);
        }

        FF_element::apply(ele, tb);
    }

    simple_timer_stop("trigon propagate");
    std::cout << "\ntrigon propagated\n";

    tb.checkout_particles();
    auto jac = tb.get_jacobian(0);

    screen << std::setprecision(16);

    for(int i=0; i<6; ++i) {
        for(int j=0; j<6; ++j) {
            screen << jac(i, j) << ", ";
        }
        screen << "\n";
    }
    screen << "\n";

    simple_timer_print(screen);
    return;
}

int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    run();

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

