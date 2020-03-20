

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/split_operator_stepper_elements.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/bunch/populate.h"
//#include "synergia/bunch/diagnostics_track.h"
//#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/collective/space_charge_2d_open_hockney.h"


int run()
{
    Logger screen(0, LoggerV::INFO_TURN);

#if 0
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
#endif

    MadX_reader reader;
    auto lattice = reader.get_lattice("oqo", "fodo.madx");
    auto const & ref = lattice.get_reference_particle();

    lattice.print(screen);

    // space charge
    Space_charge_2d_open_hockney_options sc_ops(32, 32, 128);
    sc_ops.comm_group_size = 1;

    // stepper
    Split_operator_stepper_elements stepper(sc_ops, 1);

    // Propagator
    Propagator propagator(lattice, stepper);
    //propagator.print_steps(screen);

    // bunch simulator
    auto sim = Bunch_simulator::create_single_bunch_simulator(
            lattice.get_reference_particle(), 1024 * 1024 * 1, 1e13,
            Commxx() );

    screen << "reference momentum = " << ref.get_momentum() << " GeV\n";

    // populate particle data
    auto & bunch = sim.get_bunch();
    auto local_num = bunch.get_local_num();
    auto hparts = bunch.get_host_particles();

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
    //populate_6d(dist, bunch, means, covariances);

    bunch.checkout_particles();

    double sum = 0;
    for(int p=0; p<bunch.get_local_num(); ++p)
        for(int i=0; i<6; ++i) sum += hparts(p, i);

    screen(LoggerV::DEBUG) << std::setprecision(8)
               << "\n\npopulated sum = " << sum << "\n";

    for (int p=0; p<4; ++p) bunch.print_particle(p, screen);
    screen << "\n";

#if 0
    // diagnostics
    Diagnostics_track diag_track(2, "part_2_track.h5");
    sim.reg_diag_per_turn("track_2", diag_track);

    Diagnostics_bulk_track diag_bulk_track(6, 0, "bulk_track.h5");
    sim.reg_diag_per_turn("bulk_track", diag_bulk_track);
#endif

    // propagate
    propagator.propagate(sim, screen, 4);

    // print particles after propagate
    bunch.checkout_particles();

    sum = 0;
    for(int p=0; p<bunch.get_local_num(); ++p)
        for(int i=0; i<6; ++i) sum += hparts(p, i);

    screen(LoggerV::DEBUG) << std::setprecision(8)
               << "\n\npropagated sum = " << sum << "\n";

    for (int p=0; p<4; ++p) bunch.print_particle(p, screen);

    double g_sum = 0;
    MPI_Reduce(&sum, &g_sum, 1, MPI_DOUBLE, MPI_SUM, 0, bunch.get_comm());

    screen(LoggerV::DEBUG) << std::setprecision(8)
               << "\n\npropagated sum (reduced) = " << g_sum << "\n";

    screen << "\n";

#if 0
    auto sub = bunch.get_particles_in_range(2, 3);
    for(int i=0; i<3; ++i)
    {
        screen(LoggerV::DEBUG) << std::setprecision(8)
            << sub(i, 0) << ", "
            << sub(i, 1) << ", "
            << sub(i, 2) << ", "
            << sub(i, 3) << ", "
            << sub(i, 4) << ", "
            << sub(i, 5) << "\n";
    }
#endif


    return 0;
}

template<typename T>
void print_layout_props(T const& a)
{
    std::cout 
        << "dim: " << a.extent(0) << ", " << a.extent(1)
        << ", stride: " << a.stride(0) << ", " << a.stride(1)
        << ", size = " << a.size()
        << ", cap = " << a.span()
        << "\n";
}

void layout_test()
{
    std::cout << "\nkarray1d_dev allow padding\n";
    for (int x = 10; x<20; ++x)
    {
        karray1d_dev b(Kokkos::view_alloc("b", Kokkos::AllowPadding), x);
        print_layout_props(b);
    }

    std::cout << "\nkarray1d_dev no padding\n";
    for (int x = 10; x<20; ++x)
    {
        karray1d_dev b(Kokkos::view_alloc("b"), x);
        print_layout_props(b);
    }


    std::cout << "\nkarray2d_dev allow padding\n";
    for (int x = 10; x<20; ++x)
    {
        karray2d_dev a(Kokkos::view_alloc("a", Kokkos::AllowPadding), x, 7);
        print_layout_props(a);
    }

    std::cout << "\nkarray2d_dev no padding\n";
    for (int x = 10; x<20; ++x)
    {
        karray2d_dev a(Kokkos::view_alloc("a"), x, 7);
        print_layout_props(a);
    }

    std::cout << "\nkarray2d_row_dev allow padding\n";
    for (int x = 10; x<20; ++x)
    {
        karray2d_row_dev a(Kokkos::view_alloc("a", Kokkos::AllowPadding), x, 7);
        print_layout_props(a);
    }

    std::cout << "\nkarray2d_row_dev no padding\n";
    for (int x = 10; x<20; ++x)
    {
        karray2d_row_dev a(Kokkos::view_alloc("b"), x, 7);
        print_layout_props(a);
    }

    std::cout << "\n";
}


int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);
    Kokkos::initialize(argc, argv);

    layout_test();
    run();

    Kokkos::finalize();
    MPI_Finalize();
    return 0;
}

