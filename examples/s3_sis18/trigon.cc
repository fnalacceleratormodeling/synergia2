
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

    // tunes and cdt
    auto res = Lattice_simulator::calculate_tune_and_cdt(lattice, 0);

    std::cout << "tune_h = " << res[0] << "\n";
    std::cout << "tune_v = " << res[1] << "\n";
    std::cout << "c_delta_t = " << res[2] << "\n\n";

    // chromaticities
    double dpp = 0.01;
    auto chroms = Lattice_simulator::get_chromaticities(lattice, dpp);

    std::cout << "slip_factor = " << chroms.slip_factor << "\n";
    std::cout << "slip_factor_p = " << chroms.slip_factor_prime << "\n";
    std::cout << "momentum_compaction = " << chroms.momentum_compaction << "\n";
    std::cout << "h_chrom = " << chroms.horizontal_chromaticity << "\n";
    std::cout << "h_chrom_p = " << chroms.horizontal_chromaticity_prime << "\n";
    std::cout << "v_chrom = " << chroms.vertical_chromaticity << "\n";
    std::cout << "v_chrom_p = " << chroms.vertical_chromaticity_prime << "\n\n";

    // closed orbit
    auto probe = Lattice_simulator::calculate_closed_orbit(lattice);

    std::cout << "closed orbit probe: ";
    for(int i=0; i<6; ++i) std::cout << probe[i] << ", ";
    std::cout << "\n";

    // trigon bunch
    using trigon_t = Trigon<double, 1, 6>;

    Commxx comm;
    bunch_t<trigon_t> tb(ref, comm.size(), comm);

    auto parts = tb.get_host_particles();

    // init value
    for(int i=0; i<6; ++i) 
        parts(0, i).set(probe[i], i);

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

