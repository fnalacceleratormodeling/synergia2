#include <iostream>
#include <iomanip>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/utils/serialization.h"
#include "synergia/simulation/lattice_simulator.h"

#include "beamline/beamline.h"
#include "beamline/BmlPtr.h"
#include "beamline/bmlnElmnt.h"
//#include "beamline/ElmPtr.h"
#include "beamline/Particle.h"
#include "beamline/JetParticle.h"
#include <beamline/RefRegVisitor.h>
#include <mxyzptlk/TJetVector.h>

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    const int order = 3;
    //const int order = 1;

    JetParticle::createStandardEnvironments( order );

    MadX_reader madx_reader;
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("fodo", "channel.seq"));

    std::cout << "Read lattice " << lattice_sptr->get_name() << std::endl;
    std::cout << "  lattice length: " << lattice_sptr->get_length() << std::endl;

    const Reference_particle reference_particle(lattice_sptr->get_reference_particle());

    double beam_momentum = reference_particle.get_momentum();
    double beam_energy = reference_particle.get_total_energy();

    std::cout << "Beam particle momentum: " << beam_momentum << std::endl;;
    std::cout << "Beam particle energy: " << beam_energy << std::endl;

    Lattice_simulator lattice_simulator(lattice_sptr, order);
    Chef_lattice_sptr chef_lattice_sptr(lattice_simulator.get_chef_lattice_sptr());
    BmlPtr bmlptr(chef_lattice_sptr->get_beamline_sptr());
    double brho = chef_lattice_sptr->get_brho();
    std::cout << "brho: " << brho << std::endl;


    // everything from here is just CHEF

    // Count elements
    // --------------
    std::cout << "CHEF beamline contains " << bmlptr->countHowMany() << " top level elements" << std::endl;
    std::cout << "CHEF beamline contains " << bmlptr->countHowManyDeeply() << " elements" << std::endl;

    Proton my_proton(beam_energy);
    JetProton my_jetproton( my_proton );

    std::cout << "proton initial state: " << my_proton.State() << std::endl;
    bmlptr->propagate(my_proton);
    std::cout << "proton final state: " << my_proton.State() << std::endl;
    std::cout << "=================\n" << std::endl;
    my_proton.setStateToZero();

    bmlptr->propagate(my_jetproton);

    std::cout << "map" << std::endl;
    std::cout << my_jetproton.State().jacobian() << std::endl;

    std::cout << std::endl;
    // my_jetproton.State().printCoeffs();
    TJetVector<double> mapping(my_jetproton.State());
    //mapping.printCoeffs();
    
    std::cout << "mapping dimension: " << mapping.Dim() << std::endl;
    std::cout << "mapping weight (order): " << mapping.Weight() << std::endl;

    const char *coord_names[] = {"x", "y", "cdt", "px", "py", "dpop"};

    EnvPtr<double> env = mapping.Env();
    for (int i=0; i<6; ++i) {
        std::cout << "\ncoordinate " << coord_names[i] << std::endl;
        int cnt = 0;
        for (TJet<double>::const_iterator jet_it = mapping(i).begin(); jet_it != mapping(i).end(); ++jet_it) {
            if (jet_it->coefficient() == 0.0) {
                continue;
            }
            if (cnt) {
                    std::cout << "\n";
            }
            std::cout << std::setprecision(16) << jet_it->coefficient() << " ";
            for (int j=0; j<jet_it->exponents(env).Dim(); ++j) {
                if (jet_it->exponents(env)(j)) {
                    std::cout << coord_names[j] << "^" << jet_it->exponents(env)(j) << " ";
                    //std::cout << jet_it->exponents(env);
                }
            }
            ++cnt;
        }
        std::cout << std::endl;
    }


    #if 0
    my_proton.set_x( 0.001 );
    my_proton.set_y( 0.001 );
    std::cout << "proton initial state: " << my_proton.State() << std::endl;

    for( int i = 0; i < 100; ++i ) 
    {
      bmlptr->propagate(my_proton);
      std::cout << "proton final state: " << my_proton.State() << std::endl;
    }
    my_proton.setStateToZero();
    #endif
    

    #if 0
    std::cout << "Propagation of proton through lattice" << std::endl;
    my_proton.setStateToZero();
    for ( beamline::deep_iterator it = bmlptr->deep_begin(); it != bmlptr->deep_end(); ++it ) 
    {
      std::cout << (*it)->Type() << "  " << (*it)->Name() << std::endl;
      std::cout << "In  state = " << my_proton.State() << std::endl;
      (*it)->propagate( my_proton );
      std::cout << "Out state = " << my_proton.State() << std::endl;
    }
    my_proton.setStateToZero();
    #endif

    #if 0
    std::cout << "Propagation of proton through lattice" << std::endl;
    std::cout << "Excluding thinpoles" << std::endl;
    my_proton.setStateToZero();
    for ( beamline::deep_iterator it = bmlptr->deep_begin(); it != bmlptr->deep_end(); ++it ) 
    {
      if( std::string((*it)->Type()) != std::string("ThinPole") ) {
        std::cout << (*it)->Type() << "  " << (*it)->Name() << std::endl;
        std::cout << "In  state = " << my_proton.State() << std::endl;
        (*it)->propagate( my_proton );
        std::cout << "Out state = " << my_proton.State() << std::endl;
      }
      else {
        std::cout << "ThinPole " 
                  << boost::dynamic_pointer_cast<ThinPole>(*it)->getPole() 
                  << "  " << (*it)->Name() 
                  << " strength = " << (*it)->Strength() << std::endl;
      }
    }
    my_proton.setStateToZero();

    std::cout << "Propagation of jetproton through lattice" << std::endl;
    std::cout << "Excluding thinpoles" << std::endl;
    for ( beamline::deep_iterator it = bmlptr->deep_begin(); it != bmlptr->deep_end(); ++it ) 
    {
      if( std::string((*it)->Type()) != std::string("ThinPole") ) {
        std::cout << (*it)->Type() << "  " << (*it)->Name() << std::endl;

        std::cout << "In  state = ( ";
        for( int i = 0; i < 5; ++i ) {
          std::cout << my_jetproton.State()[i].standardPart() << ", ";
        }
        std::cout << my_jetproton.State()[5].standardPart() << " )";
        std::cout << std::endl;

        (*it)->propagate( my_jetproton );

        std::cout << "Out state = ( ";
        for( int i = 0; i < 5; ++i ) {
          std::cout << my_jetproton.State()[i].standardPart() << ", ";
        }
        std::cout << my_jetproton.State()[5].standardPart() << " )";
        std::cout << std::endl;
      }
      else {
        std::cout << "ThinPole " 
                  << boost::dynamic_pointer_cast<ThinPole>(*it)->getPole() 
                  << "  " << (*it)->Name() 
                  << " strength = " << (*it)->Strength() << std::endl;
      }
    }
    #endif



    #if 0
    bool toggled = false;
    std::cout << "\n\nHere we go!\n" << std::endl;
    // bmlptr->propagate( my_jetproton );
    for ( beamline::deep_iterator it = bmlptr->deep_begin(); it != bmlptr->deep_end(); ++it ) 
    {
      if(toggled) {
        std::cout << "Toggled element" << std::endl;
        std::cout << (*it)->Type() << "  " << (*it)->Name() << std::endl;
        std::cout << "length = " << (*it)->Length() << "; strength = " << (*it)->Strength() << std::endl;

        std::cout << "In state: " << endl;
        my_jetproton.State().printCoeffs();

        (*it)->propagate( my_jetproton );

        std::cout << "Out state: " << endl;
        my_jetproton.State().printCoeffs();
      }
      else {
        if( (*it)->Name() == std::string("PR.BHU000028.mpdi1_12pole" ) ) {
          toggled = true;
        }
        std::cout << (*it)->Type() << "  " << (*it)->Name() << std::endl;
        (*it)->propagate( my_jetproton );
      }
    }
    #endif
}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    run();
    MPI_Finalize();
    return 0;
}

