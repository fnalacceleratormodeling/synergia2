
#include "synergia/foundation/physical_constants.h"

#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"

#include "synergia/simulation/split_operator_stepper.h"
//#include "synergia/simulation/independent_stepper_elements.h"

#include "synergia/bunch/populate.h"
#include "synergia/bunch/core_diagnostics.h"

#include "synergia/lattice/madx_reader.h"
#include "synergia/utils/lsexpr.h"

#include "synergia/collective/space_charge_2d_open_hockney.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"
//#include "synergia/collective/dummy_collective_operator.h"

#include "synergia/utils/simple_timer.h"


void print_bunch_statistics(Bunch_sptr bunch, int rank)
{
    MArray1d mean(Core_diagnostics::calculate_mean(*bunch));
    MArray1d std(Core_diagnostics::calculate_std(*bunch, mean));

    if (rank !=0 ) return;

    std::cerr
         << std::setprecision(16)
         << std::showpos
         << std::ios::scientific;

    std::cerr << std::endl;

    for(int i=0; i<6; ++i)
        std::cerr << i << ": " << mean[i] << " \t " << std[i] << "\n";

    std::cerr << std::endl;

    auto parts = bunch->get_local_particles();

    std::cerr << "Particles 1-4" << std::endl;
    for(int p=0; p<4; ++p) 
    {
        std::cerr << p << ":";
        for(int i=0; i<6; ++i)
            std::cerr << "  " << parts[p][i];
        std::cerr << "\n";
    }

    std::cerr << "\n";
}


void run_and_save(std::string & prop_str, std::string & sim_str)
{
    Commxx_sptr comm(new Commxx());
    const int myrank = comm->get_rank();
    const int mpisize = comm->get_size();

    std::vector<int> grid_shape = {64, 64, 64};

    static std::string lattice_file = "sis18-6.lsx";
    auto lsexpr = read_lsexpr_file(lattice_file);
    Lattice_sptr lattice(new Lattice(lsexpr));

    lattice->set_all_string_attribute("extractor_type", "libff");

    if (myrank == 0) {
        std::cout << "read lattice from file " << lattice_file << std::endl;
        std::cout << "    length: " << lattice->get_length() << std::endl;
        std::cout << "    # elements: " << lattice->get_elements().size() << std::endl;
    }

    Lattice_simulator lsim(lattice, 2);

    if (myrank == 0) {
        std::cout << "created lattice simulator" << std::endl;
    }

    // tune the lattice
    //Lattice_simulator::tune_circular_lattice(lattice);
    //lsim.tune_circular_lattice();

    // get the reference particle
    auto const & ref = lattice->get_reference_particle();

    if (myrank == 0) {
        std::cout << "reference momentum = " 
                  << ref.get_momentum() << " GeV\n";
    }

    // space charge
    //Commxx_divider_sptr commxx_divider_sptr(new Commxx_divider(32, false));
    //Space_charge_2d_open_hockney_sptr sc(new Space_charge_2d_open_hockney(commxx_divider_sptr, {64, 64, 64}));
    Space_charge_3d_open_hockney_sptr sc(new Space_charge_3d_open_hockney({64, 64, 64}));
    //Dummy_collective_operator_sptr sc(new Dummy_collective_operator("dummy"));
    if (myrank == 0) {
        std::cout << "created space charge operator: " << sc->get_type() << ":" << sc->get_name() << std::endl;
    }

    // stepper
    Split_operator_stepper_sptr stepper(
            new Split_operator_stepper(lsim, sc, 71));

    if (myrank == 0) {
        std::cout << "created stepper" << std::endl;
    }

    // Propagator
    Propagator propagator(stepper);

    // bunch simulator
    Bunch_sptr bunch(new Bunch(ref, 4194394, 2.94e10, comm));
    bunch->set_bucket_index(0);

    const std::string bunchdatafile = "turn_particles_0000.h5";
    if (myrank == 0) {
        std::cout << "Reading bunch data from file: " << bunchdatafile << "...";
    }

    // or read from file
    bunch->read_file(bunchdatafile);
    if (myrank == 0) {
        std::cout << " read " << bunch->get_total_num() << " particles" << std::endl;
    }

    // bunch statistics
    if (myrank == 0) {
        std::cout << "Initial bunch statistics" << std::endl;
    }

    print_bunch_statistics(bunch, myrank);

    if (myrank == 0) {
        std::cout << std::endl;
    }

    // bunch simulator
    Bunch_simulator bsim(bunch);

    // diagnostics
    auto diag_loss = Diagnostics_loss_sptr(new Diagnostics_loss("loss.h5", "aperture"));
    diag_loss->set_bunch_sptr(bsim.get_bunch_sptr());
    lsim.get_lattice_sptr()->add_loss_diagnostics(diag_loss);

    // propagate
    propagator.propagate(bsim, 1, 1, 2);

    // bunch statistics
    if (myrank == 0) {
        std::cout << std::endl;
        std::cout << "Final bunch statistics" << std::endl;
    }

    print_bunch_statistics(bunch, myrank);

    return;
}

int main(int argc, char ** argv)
{
    MPI_Init(&argc, &argv);

    std::string prop_str, sim_str;
    run_and_save(prop_str, sim_str);

    MPI_Finalize();
    return 0;
}

