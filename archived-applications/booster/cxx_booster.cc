#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/split_operator_stepper_choice.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    std::cout<<" run start"<<std::endl;
    const double num_macroparticles=100000;  
    const int map_order = 1;
    const int num_steps=144;
    const bool bpms=true;
    
    Lattice_sptr lattice_sptr(new Lattice());
    std::cout<<"before lattice load "<<std::endl;
    try {
        xml_load(*lattice_sptr, "booster_lattice.xml");
    }
    catch (std::runtime_error) {
        std::cerr << "cxx_booster: failed to find booster_lattice.xml\n";
        std::cerr << "Run booster_xml.py to generate booster_lattice.xml\n";
        exit(1);
    }

    std::cout<<"lattice loaded "<<std::endl;

   

// ***************************************
/*
  double ref[1]={1.};
  double scale[1]={1.e-3};
  double * pref=ref;
  double * pscale=scale;
  int maxweight=1;
  int numvar=1;
  int spacedim=1;

  EnvPtr<double> evp=
  TJetEnvironment<double>::makeJetEnvironment(maxweight, numvar, spacedim, pref,pscale);*/
// ***************************************

   Chef_lattice_sptr  chef_lattice_sptr(new Chef_lattice(lattice_sptr));
   std::cout<<"chef lattice constructed "<<std::endl;
   BmlPtr beamline_sptr(chef_lattice_sptr->get_beamline_sptr()->Clone());
   std::cout<<"beamline_sptr created "<<std::endl;
   Proton probe;
   double momentum(lattice_sptr->get_reference_particle().get_momentum());
   probe.SetReferenceMomentum(momentum);
   probe.setStateToZero();

   if (Jet__environment::getLastEnv() == 0) {
           JetParticle::createStandardEnvironments(map_order);
   }

   BeamlineContext probecontext(probe, beamline_sptr); 
   probecontext.handleAsRing();
   std::cout<<"beamlinecontext created "<<std::endl;


  
    double tune=probecontext.getHorizontalEigenTune();
    std::cout<<" tune="<<tune<<std::endl;
// ***************************************




//    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
//      std::cout<<"lattice simulator passed "<<std::endl;

//    double chef_frac_tunex=lattice_simulator.get_horizontal_tune();
//    double chef_fracnedit _tuney=lattice_simulator.get_vertical_tune();
//       std::cout<<" chef_frac_tunex="<<chef_frac_tunex<<std::endl;
//        std::cout<<" chef_frac_tuney="<<chef_frac_tuney<<std::endl;
//    double chef_eigen_tuqnex=lattice_simulator.get_horizontal_tune(true);
//    double chef_eigen_tuney=lattice_simulator.get_vertical_tune(true);
//        std::cout<<" chef_eigen_tunex="<<chef_eigen_tunex<<std::endl;
//        std::cout<<" chef_eigen_tuney="<<chef_eigen_tuney<<std::endl;
//    double horizontal_chromaticity=lattice_simulator.get_horizontal_chromaticity();
//       std::cout<<"horizontal chromaticity done"<<std::endl;
//    double vertical_chromaticity=lattice_simulator.get_vertical_chromaticity();
//    double momentum_compaction=lattice_simulator.get_momentum_compaction();
//    double slip_factor=lattice_simulator.get_slip_factor();
//        std::cout<<" horizontal_chromaticity="<<horizontal_chromaticity<<std::endl;
//        std::cout<<" vertical_chromaticity="<<vertical_chromaticity<<std::endl;
//        std::cout<<" momentum_compaction="<<momentum_compaction<<std::endl;
//        std::cout<<" slip_factor="<<slip_factor<<std::endl;
/*
    Lattice_elements  horizontal_correctors;
    Lattice_elements  vertical_correctors;

   for (Lattice_elements::const_iterator latt_it =
        lattice_simulator.get_lattice().get_elements().begin();
        latt_it !=lattice_simulator.get_lattice().get_elements().end(); ++latt_it){
      
        if ((*latt_it)->get_name().find("ssxs")==0){
             //  std::cout<<" horizontal sextupole correctors: "<<std::endl;
             //  std::cout<<" element_name: "<<(*latt_it)->get_name()<<std::endl;
               horizontal_correctors.push_back(*latt_it);
        }
       
       if ((*latt_it)->get_name().find("ssxl")==0){
         // std::cout<<" vertical sextupole correctors :"<<std::endl; 
         //     std::cout<<" element_name: "<<(*latt_it)->get_name()<<std::endl;
              vertical_correctors.push_back(*latt_it);
      }

    }
   std::cout<<"correctors done"<<std::endl;

   double const chromH=-10.;
   double const chromV=-4.;
   double const tolerance=1.e-3;
   int const max_steps=10;

   lattice_simulator.adjust_chromaticities(chromH, chromV, 
        horizontal_correctors, vertical_correctors, tolerance, max_steps);
    
        std::cout<<"adjust_chromaticities done"<<std::endl;

//    horizontal_chromaticity=lattice_simulator.get_horizontal_chromaticity();
//    vertical_chromaticity=lattice_simulator.get_vertical_chromaticity();
//    momentum_compaction=lattice_simulator.get_momentum_compaction();
//    slip_factor=lattice_simulator.get_slip_factor();
//        std::cout<<" horizontal_chromaticity="<<horizontal_chromaticity<<std::endl;
//        std::cout<<" vertical_chromaticity="<<vertical_chromaticity<<std::endl;
//        std::cout<<" momentum_compaction="<<momentum_compaction<<std::endl;
//        std::cout<<" slip_factor="<<slip_factor<<std::endl;

Stepper_sptr stepper_sptr;

if (bpms){
    List_choice_map list_choice_map;

    const int steps_per_bpm=1;
    Dummy_collective_operator_sptr bpm_measure_sptr(new Dummy_collective_operator("bmp_measure"));
    std::cout<<" dummy bmp_measure "<<std::endl;

    Collective_operators operators_bpm;
    operators_bpm.push_back(bpm_measure_sptr);
    Kicks  kicks_bpm(operators_bpm,steps_per_bpm);

    for (Lattice_elements::const_iterator latt_it =
        lattice_sptr->get_elements().begin();
        latt_it !=lattice_sptr->get_elements().end(); ++latt_it){
      
        if ( !((*latt_it)->get_name().find("obpmold")==0) &&
                ((*latt_it)->get_name().find("obpm")==0)  ){             
              // std::cout<<" element_name: "<<(*latt_it)->get_name()<<std::endl;
            list_choice_map[(*latt_it)->get_name()]=kicks_bpm;
        }
    }

    const int num_steps_else=96;
   Dummy_collective_operator_sptr noop_sptr(new Dummy_collective_operator("stub"));

    Collective_operators  operators_else;
    operators_else.push_back(noop_sptr);
    Kicks kicks_else(operators_else,num_steps_else);
    list_choice_map["else"]=kicks_else;
    std::cout<<" list_choice_map done"<<std::endl;
    stepper_sptr=Stepper_sptr(new Split_operator_stepper_choice(lattice_sptr, map_order, list_choice_map, true ));
    std::cout<<" split choice stepper done"<<std::endl;
}
else{
  stepper_sptr=Stepper_sptr(new  Independent_stepper(lattice_sptr, map_order, num_steps));
  std::cout<<" independent stepper done"<<std::endl;
}

    double horizontal_chromaticity=stepper_sptr->get_lattice_simulator().get_horizontal_chromaticity();
    double  vertical_chromaticity=stepper_sptr->get_lattice_simulator().get_vertical_chromaticity();
    double  momentum_compaction=stepper_sptr->get_lattice_simulator().get_momentum_compaction();
    double  slip_factor=stepper_sptr->get_lattice_simulator().get_slip_factor();
        std::cout<<" horizontal_chromaticity="<<horizontal_chromaticity<<std::endl;
        std::cout<<" vertical_chromaticity="<<vertical_chromaticity<<std::endl;
        std::cout<<" momentum_compaction="<<momentum_compaction<<std::endl;
        std::cout<<" slip_factor="<<slip_factor<<std::endl;

    std::cout<<" run done"<<std::endl;   */ 
}
int
main(int argc, char **argv)
{
  // MPI_Init(&argc, &argv);
  run();
  // MPI_Finalize();
    return 0;
}
