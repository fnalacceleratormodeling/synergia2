
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


void print_bunch_statistics(Bunch_sptr bunch)
{
    MArray1d mean(Core_diagnostics::calculate_mean(*bunch));
    MArray1d std(Core_diagnostics::calculate_std(*bunch, mean));

    std::cerr
         << std::setprecision(16)
         << std::showpos
         << std::ios::scientific;

    for(int i=0; i<6; ++i)
        std::cerr << mean[i] << ", " << std[i] << "\n";

    auto parts = bunch->get_local_particles();

    for(int p=0; p<4; ++p) 
    {
        for(int i=0; i<6; ++i)
            std::cerr << parts[p][i] << ", ";
        std::cerr << "\n";
    }

    std::cerr << "\n";
}


void run_and_save(std::string & prop_str, std::string & sim_str)
{
    std::vector<int> grid_shape = {64, 64, 64};

    auto lsexpr = read_lsexpr_file("sis18-6.lsx");
    Lattice_sptr lattice(new Lattice(lsexpr));

    lattice->set_all_string_attribute("extractor_type", "libff");

    Lattice_simulator lsim(lattice, 2);

    // tune the lattice
    //Lattice_simulator::tune_circular_lattice(lattice);
    //lsim.tune_circular_lattice();

    // get the reference particle
    auto const & ref = lattice->get_reference_particle();

    std::cout << "reference momentum = " 
        << ref.get_momentum() << " GeV\n";

    // space charge
    //Commxx_divider_sptr commxx_divider_sptr(new Commxx_divider(32, false));
    //Space_charge_2d_open_hockney_sptr sc(new Space_charge_2d_open_hockney(commxx_divider_sptr, {64, 64, 64}));
    Space_charge_3d_open_hockney_sptr sc(new Space_charge_3d_open_hockney({64, 64, 64}));
    //Dummy_collective_operator_sptr sc(new Dummy_collective_operator("dummy"));

    // stepper
    Split_operator_stepper_sptr stepper(
            new Split_operator_stepper(lsim, sc, 71));

    // Propagator
    Propagator propagator(stepper);

    // bunch simulator
    Commxx_sptr comm(new Commxx());
    Bunch_sptr bunch(new Bunch(ref, 4194394, 2.94e10, comm));
    bunch->set_bucket_index(0);

    // or read from file
    bunch->read_file("turn_particles_0000_4M.h5");

    // bunch statistics
    print_bunch_statistics(bunch);

    // bunch simulator
    Bunch_simulator bsim(bunch);

    // diagnostics
    auto diag_loss = Diagnostics_loss_sptr(new Diagnostics_loss("loss.h5", "aperture"));
    diag_loss->set_bunch_sptr(bsim.get_bunch_sptr());
    lsim.get_lattice_sptr()->add_loss_diagnostics(diag_loss);

    // propagate
    propagator.propagate(bsim, 1, 1, 2);

    // bunch statistics
    print_bunch_statistics(bunch);

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

