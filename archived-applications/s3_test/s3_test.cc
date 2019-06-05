

#include "synergia/simulation/propagator.h"


int run()
{
    Logger screen(0, LoggerV::DEBUG);

    Lattice lattice("test");

    Lattice_element e1("drift", "d1");
    e1.set_double_attribute("l", 1.0);
    lattice.append(e1);

    Lattice_element e2("drift", "d2");
    e2.set_double_attribute("l", 0.8);
    lattice.append(e2);

    lattice.print();

    Propagator propagator(lattice);
    propagator.print_steps(screen);

    return 0;
}


int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);

    run();

    MPI_Finalize();
    return 0;
}

