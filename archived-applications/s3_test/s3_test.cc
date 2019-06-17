

#include "synergia/simulation/propagator.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/populate.h"



int run()
{
    Logger screen(0, LoggerV::DEBUG);

    // Lattice
    Lattice lattice("test");

    Reference_particle ref(1, pconstants::mp, 3.0);
    lattice.set_reference_particle(ref);

    Lattice_element e1("drift", "d1");
    e1.set_double_attribute("l", 1.0);
    lattice.append(e1);

    Lattice_element e2("quadrupole", "q1");
    e2.set_double_attribute("l", 1.0);
    e2.set_double_attribute("k1", 0.02);
    e2.set_double_attribute("k1s", 0.03);
    lattice.append(e2);

    Lattice_element e3("drift", "d2");
    e3.set_double_attribute("l", 0.8);
    lattice.append(e3);

    lattice.print(screen);

    // Propagator
    Propagator propagator(lattice);
    propagator.print_steps(screen);

    // bunch simulator
    auto sim = Bunch_simulator::create_single_bunch_simulator(
            lattice.get_reference_particle(), 1024, 1e13 );

    screen << ref.get_momentum() << "\n";

    // init particle data
    auto & bunch = sim.get_bunch();
    auto local_num = bunch.get_local_num();
    auto hparts = bunch.get_host_particles();

#if 0
    bunch.checkout_particles();

    for (int p=0; p<local_num; ++p)
    {
        for (int i=0; i<6; ++i)
        {
            hparts(p, i) = p*0.1 + i*0.01;
        }
    }

    bunch.checkin_particles();
#endif

    karray1d means("means", 6);
    for (int i=0; i<6; ++i) means(i) = 0.0;

    karray2d covariances("covariances", 6, 6);
    for (int i=0; i<6; ++i)
        for (int j=0; j<6; ++j)
            covariances(i, j) = 0.0;

    covariances(0,0) = 1e-6;
    covariances(1,1) = 1e-6;
    covariances(2,2) = 1e-6;
    covariances(3,3) = 1e-6;
    covariances(4,4) = 1e-6;
    covariances(5,5) = 1e-6;

    Random_distribution dist(5, Commxx());
    populate_6d(dist, bunch, means, covariances);

    bunch.checkout_particles();

    screen(LoggerV::DEBUG) << "\n\npopulated\n";
    for (int p=0; p<4; ++p) bunch.print_particle(p, screen);
    screen << "\n";

    // propagate
    sim.set_turns(0, 1);
    propagator.propagate(sim, screen);

    bunch.checkout_particles();

    screen(LoggerV::DEBUG) << "\n\npropagated\n";
    for (int p=0; p<4; ++p) bunch.print_particle(p, screen);
    screen << "\n";

    return 0;
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

