#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/commxx_per_host.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/stepper.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/collective/ecloud_from_vorpal.h"
//
// test example for e-cloud collective effect. Cloned and adapted from cxx_test.cc 
// 
 
// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run(bool do_space_charge, bool do_ecloud, const std::string &file_name_ecloud)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 32;
    grid_shape[1] = 32;
    grid_shape[2] = 256;
    const int part_per_cell = 10;
    const int num_macro_particles = grid_shape[0] * grid_shape[1]
            * grid_shape[2] * part_per_cell;
    const int seed = 4;
    const double num_real_particles = 3e14;
    const int num_steps = 8;
    const int num_turns = 30;
    const int map_order = 2;

    Lattice_sptr lattice_sptr(new Lattice());
    try {
        xml_load(*lattice_sptr, "cxx_lattice.xml");
    }
    catch (std::runtime_error) {
        std::cerr << "cxx_ecloudexample: failed to find cxx_lattice.xml\n";
        std::cerr << "Run cxx_ecloudexample.py to generate cxx_lattice.xml\n";
        exit(1);
    }
    Stepper_sptr stepper_sptr;
    
    Commxx_sptr commxx_per_host_sptr(new Commxx(true));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    std::string file_name_out("cxx_eCloud_");
    
    if (do_space_charge && (!do_ecloud)) {
      Space_charge_3d_open_hockney_sptr space_charge_sptr(
            new Space_charge_3d_open_hockney(commxx_per_host_sptr, grid_shape));
      space_charge_sptr->set_charge_density_comm(Space_charge_3d_open_hockney::charge_allreduce);
      if (num_steps < 0) {
        file_name_out += std::string("SpaceChargeMDElem3e14");
        stepper_sptr = Split_operator_stepper_elements_sptr(new Split_operator_stepper_elements(lattice_simulator, space_charge_sptr, 1));
      } else {
        file_name_out += std::string("SpaceChargeMDFixed3e14");      
        stepper_sptr = Split_operator_stepper_sptr(new Split_operator_stepper(lattice_simulator, space_charge_sptr, num_steps));
      }
    } else if ((!do_space_charge) && (do_ecloud)) {
      file_name_out += std::string("ECloudxxx");
      std::cerr << " About to instantiate the collective operator for e cloud from file " << file_name_ecloud << std::endl;
      Ecloud_from_vorpal_sptr  e_cloud_sptr(new Ecloud_from_vorpal(commxx_per_host_sptr, file_name_ecloud));
      stepper_sptr = Split_operator_stepper_elements_sptr(new Split_operator_stepper_elements(lattice_simulator, e_cloud_sptr, num_steps));
    }else if (do_space_charge && do_ecloud) {    
      file_name_out += std::string("SpaceChargeMDECloudxxx");
      Space_charge_3d_open_hockney_sptr space_charge_sptr(
            new Space_charge_3d_open_hockney(commxx_per_host_sptr, grid_shape));
      space_charge_sptr->set_charge_density_comm(Space_charge_3d_open_hockney::charge_allreduce);
      Ecloud_from_vorpal_sptr  e_cloud_sptr(new Ecloud_from_vorpal(commxx_per_host_sptr, file_name_ecloud, std::string("quadrupole")));
      Collective_operators two_ops;
      two_ops.push_back(space_charge_sptr);
      two_ops.push_back(e_cloud_sptr);       
//      stepper_sptr = Split_operator_stepper_sptr(new Split_operator_stepper(lattice_simulator, two_ops, num_steps));
// Jan 9 2013: try the split operator elements.  Works, but call the ecould operator only once per elements. 
// To be discussed!. 
//
      stepper_sptr = Split_operator_stepper_elements_sptr(new Split_operator_stepper_elements(lattice_simulator, two_ops, 1));
    } else {
       file_name_out += std::string("none");
       stepper_sptr = Stepper_sptr(new Stepper(lattice_simulator));
     }	    
     
    file_name_out += std::string(".h5");
    	    
    Propagator propagator(stepper_sptr);

    Commxx_sptr comm_sptr(new Commxx);
    Bunch_sptr bunch_sptr(
            new Bunch(lattice_sptr->get_reference_particle(),
                    num_macro_particles, num_real_particles, comm_sptr));
    Random_distribution distribution(seed, *comm_sptr);
    MArray1d means;
    xml_load(means, "cxx_means.xml");
    MArray2d covariances;
    xml_load(covariances, "cxx_covariance_matrix.xml");
    populate_6d(distribution, *bunch_sptr, means, covariances);

    Bunch_simulator bunch_simulator(bunch_sptr);
//    bunch_simulator.get_diagnostics_actions().add_per_step(
//            Diagnostics_sptr(
//                    new Diagnostics_basic("cxx_example_per_step.h5")));
    bunch_simulator.get_diagnostics_actions().add_per_turn(
            Diagnostics_sptr(
                    new Diagnostics_full2(file_name_out.c_str())));

    propagator.set_checkpoint_period(2000);
    propagator.set_final_checkpoint(true);
    double t0 = MPI_Wtime();
    const int max_turns = 0;
    const int verbosity = 1;
    propagator.propagate(bunch_simulator, num_turns, max_turns, verbosity);
    double t1 = MPI_Wtime();
    if (comm_sptr->get_rank() == 0) {
        std::cout << "propagate time = " << (t1 - t0) << std::endl;
    }
}
int
main(int argc, char **argv)
{
    bool do_space_charge = true;
    bool do_ecloud = false;    
    MPI_Init(&argc, &argv);
    std::string name=std::string("/data/lebrun/Synergia/ECloudMaps/Efield_MI2D-S2-V2f-a7.bin");    
    run(do_space_charge, do_ecloud, name);
    MPI_Finalize();
    return 0;
}
