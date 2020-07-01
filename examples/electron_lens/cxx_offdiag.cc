#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>
#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/utils/lsexpr.h"
#include <boost/archive/archive_exception.hpp>
#include "synergia/utils/multi_array_to_string.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/split_operator_stepper_elements.h"
#include "synergia/simulation/independent_stepper.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/split_operator_stepper_choice.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/fast_normal_form.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/simulation/populate_stationary.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_phase_space_density.h"
#include "synergia/simulation/diagnostics_normal_form.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/collective/space_charge_2d_open_hockney.h"
#include "synergia/collective/space_charge_rectangular.h"
#include "synergia/collective/impedance.h"

#include "cxx_offdiag_options.h"

#include "elens_actions.h"

#define DEBUG 0

//-----------------------------------------------------------------------------------

void activate_apertures(Cxx_offdiag_options const& opts, Stepper_sptr stepper_sptr, int rank, double aperture_radius)
{
    Lattice_sptr lattice_sptr(stepper_sptr->get_lattice_simulator().get_lattice_sptr());

    if (rank == 0 ) {
        std::cout << "activating aperture at " << aperture_radius << std::endl;
    }
    // activate apertures at the beginning of the cell
    int aperture_cnt = 0;
    for (Lattice_elements::iterator it=lattice_sptr->get_elements().begin();
         it!=lattice_sptr->get_elements().end(); ++it) {
        if (((*it)->get_type() == "marker") && ((*it)->get_name() == "boc")) {
            // this is the beginning of the cell.  Put a circular aperture on it
            (*it)->set_string_attribute("aperture_type", "circular");
            (*it)->set_double_attribute("circular_aperture_radius", aperture_radius);
            ++aperture_cnt;
        }
    }
    if (rank == 0) {
        std::cout << "Number of apertures set: " << aperture_cnt << std::endl;
    }
}

//-----------------------------------------------------------------------------------

void activate_elens(Cxx_offdiag_options const& opts, Stepper_sptr stepper_sptr, int rank)
{
    Lattice_sptr lattice_sptr(stepper_sptr->get_lattice_simulator().get_lattice_sptr());

    // activate electron lens
    const double elens_energy = opts.elensenergy;
    const double elens_radius = opts.elensradius;
    const double elens_length = opts.elenslength;
    double elens_current = opts.elenscurrent;
    const double elens_longrms = opts.elenslongrms;

    const int elens_divide = opts.elensdivide;

    if (rank == 0) {
        std::cout << "-----------------------------------------------------" << std::endl;
        std::cout << "Activate Electron lens" << std::endl;
        std::cout << "electron lens parameters" << std::endl;
        std::cout << "electron beam radius: " << elens_radius << std::endl;
        std::cout << "electron beam length: " << elens_length << std::endl;
        std::cout << "electron beam current: " << elens_current << std::endl;
        std::cout << "electron beam energy: " << elens_energy <<std::endl;
        std::cout << "dividing number of elens by " << elens_divide;
        std::cout << "electron pulse RMS length: ";
        if (elens_longrms <= 0.0) {
            std::cout << " UNIFORM ";
        } else {
            std::cout << elens_longrms;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }

    int elens_cnt = 0;
    for (Lattice_elements::iterator it=lattice_sptr->get_elements().begin();
         it!=lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_type() == "elens") {
            // elens needs chef_propagate
            (*it)->set_string_attribute("extractor_type", "chef_propagate");
            (*it)->set_double_attribute("radius", elens_radius);
            (*it)->set_double_attribute("gaussian", 1.0);
            (*it)->set_double_attribute("eenergy", elens_energy);
            (*it)->set_double_attribute("l", 0.0);
            (*it)->set_double_attribute("current", elens_current*elens_length);
            (*it)->set_double_attribute("longrms", elens_longrms);
            if ((elens_cnt % elens_divide) != 0) {
                // knocking out this lens
                (*it)->set_double_attribute("current", 0.0);
            }
            if (rank == 0) {
                std::cout << (*it)->as_string() << std::endl;
            }
            ++elens_cnt;
        }
    }
}

//-----------------------------------------------------------------------------------

Lattice_sptr setup_lattice(Cxx_offdiag_options const& opts, int rank)
{

    Lattice_sptr lattice_sptr(new Lattice(read_lsexpr_file(opts.lattice_file)));
    if (rank==0) {
        std::cout<<"Read lattice from file "<<opts.lattice_file<< ", " << lattice_sptr->get_elements().size() << std::endl;
        std::cout << "    lattice_length: " << lattice_sptr->get_length() << std::endl;
    }

   // determine type of propagation
    std::string propagate_type;
    if (opts.transport == "maps") {
        propagate_type = "chef_map";
    } else if (opts.transport == "chef") {
        propagate_type = "chef_propagate";
    } else if (opts.transport == "libff") {
        propagate_type = "libff";
    } else {
        throw std::runtime_error("unknown tranport type: "+opts.transport+", allowable values: maps|chef|libff");
    }

    if (rank == 0) {
        std::cout << "Selected propagation type: " << opts.transport << std::endl;
    }

    for (Lattice_elements::const_iterator it =
             lattice_sptr->get_elements().begin();
         it != lattice_sptr->get_elements().end(); ++it) {

        std::string element_name=(*it)->get_name();
        std::string element_type=(*it)->get_type();

        (*it)->set_string_attribute("extractor_type", propagate_type);

        // add forced diagnostics to cell beginning
        //~ if ((element_type == "elens") &&
            //~ (element_name == "lens1")) {
        if ((element_type == "marker") &&
            (element_name == "boc")) {
            (*it)->set_string_attribute("force_diagnostics", "true");
            if (false && (rank == 0)) {
                std::cout << "lens1: " << (*it)->as_string() << std::endl;;
            }
        }

        if ((element_type == "multipole") &&
            (element_name == "s")) {
            (*it)->set_double_attribute("k2l", opts.k2l);
            // the sextupole has to be propagated with chef
            (*it)->set_string_attribute("extractor_type", "chef_propagate");
        }
                      
        std::string rf_name(element_name.begin(),element_name.begin()+1);
        if  (rf_name=="r"){
            //std::cout<<" elemnt name= "<<(*it)->get_name()<<std::endl;
            (*it)->set_double_attribute("volt", opts.rf_voltage);
            (*it)->set_double_attribute("lag", 0.);
            (*it)->set_double_attribute("harmon", opts.harmon);
            //  (*it)->print();
        }
        //                 std::string bpm_name(element_name.begin(),element_name.begin()+3);
        //                if (bpm_name=="bpm"){
        //                     (*it)->print();
        //                }

        /*
          if (opts.chef_propagate) {
          (*it)->set_string_attribute("extractor_type", "chef_propagate");
          }
          else if (opts.chef_map) {
          (*it)->set_string_attribute("extractor_type", "chef_map");
          }

          if (opts.if_aperture) {
          std::string name1(element_name.begin(),element_name.begin()+1);
          if ((element_type=="sbend") && (name1=="m")) {
          //  (*it)->print();
          //  std::cout<<" elem name1="<<name1<<" type="<<element_type<<std::endl;
          (*it)->set_string_attribute("aperture_type","rectangular");
          (*it)->set_double_attribute("rectangular_aperture_width", 2*opts.aperture_bending);
          (*it)->set_double_attribute("rectangular_aperture_height", 2.*opts.aperture_bending);

          }
          else{
          (*it)->set_string_attribute("aperture_type","circular");
          (*it)->set_double_attribute("circular_aperture_radius", opts.aperture_straight);
          }
          (*it)->set_string_attribute("aperture_loss","aperture_loss_file");
          }*/

        //   (*it)->print();

    }

    return lattice_sptr;
}

//-----------------------------------------------------------------------------------

void print_lattice(Lattice_sptr lattice_sptr, std::string filename)
{
    std::ofstream of(filename.c_str(), std::ofstream::out);
    of << lattice_sptr->as_string() << std::endl;
    of.close();
}

//-----------------------------------------------------------------------------------

void print_beamline(Stepper_sptr stepper_sptr, std::string filename)
{
    std::ofstream of(filename.c_str(), std::ofstream::out);
    BmlPtr chef_beamlineptr(stepper_sptr->get_lattice_simulator().get_chef_lattice_sptr()->get_sliced_beamline_sptr());
    of << chef_beamline_as_string(chef_beamlineptr) << std::endl;;
    of.close();
}

//-----------------------------------------------------------------------------------

void write_lsx_lattice(Lattice_sptr lattice_sptr, std::string filename)
{
    write_lsexpr_file(lattice_sptr->as_lsexpr(),filename);
}

//-----------------------------------------------------------------------------------

void set_tunes(int rank, Stepper_sptr stepper_sptr, double nu_x, double nu_y)
{
    // need to update lattice simulator because adjust tunes is
    // going to use the values in the chef lattice so we need
    // to reflect the new values
    stepper_sptr->get_lattice_simulator().update();
    // restore the tunes to their original value
    // first collect the tuning elements
    Lattice_elements focussing, defocussing;
    Lattice_sptr lattice_sptr(stepper_sptr->get_lattice_simulator().get_lattice_sptr());
    for (Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
         it!=lattice_sptr->get_elements().end(); ++it) {
        if ((*it)->get_type() == "quadrupole") {
            double k1 = (*it)->get_double_attribute("k1");
            if (k1 > 0.0) {
                focussing.push_back(*it);
            } else if (k1 < 0.0) {
                defocussing.push_back(*it);
            }
        }
    }
    if (rank == 0) {
        std::cout << "Found " << focussing.size() << " focussing elements" << std::endl;
        std::cout << "Found " << defocussing.size() << " defocussing elements" << std::endl;
    }
    stepper_sptr->get_lattice_simulator().adjust_tunes(nu_x, nu_y, focussing, defocussing, 1.0e-6, 1);

}

//-----------------------------------------------------------------------------------

void activate_error_element(Cxx_offdiag_options const& opts, int rank, Stepper_sptr stepper_sptr, double nu_x, double nu_y)
{
    Lattice_sptr lattice_sptr(stepper_sptr->get_lattice_simulator().get_lattice_sptr());

    // find element that gets lattice errors if any
    Lattice_elements::iterator errelement_it = lattice_sptr->get_elements().end();
    int errelement = opts.errelement;
    if ((opts.errsize != 0.0) && (errelement >= 0)) {
        int ecnt = 0;
        for (Lattice_elements::iterator it=lattice_sptr->get_elements().begin();
             it!=lattice_sptr->get_elements().end(); ++it) {
            if (ecnt == errelement) {
                errelement_it = it;
                // make sure it has a k1 attribute
                if (! (*errelement_it)->has_double_attribute("k1")) {
                    throw std::runtime_error("error element does not have a k1 attribute: " + (*errelement_it)->as_string());
                }
                break;
            } else {
                ++ecnt;
            }
        }
        if (errelement_it == lattice_sptr->get_elements().end()) {
            throw std::runtime_error("error element requested, but it could not be found.");
        } else { // error element found
            if (rank == 0) {
                std::cout << "==========================================" << std::endl;
                std::cout << "Activating error element" << (*errelement_it)->get_name() << "    error size: " << opts.errsize  << std::endl;
            }
            double k1 = (*errelement_it)->get_double_attribute("k1");
            double newk1 = k1 * (1.0 + opts.errsize);
            (*errelement_it)->set_double_attribute("k1", newk1);
            if (rank == 0) {
                std::cout << "    Old k1: " << k1 << ", new k1: " << newk1 << ": " << (*errelement_it)->as_string() << std::endl;
            }

            set_tunes(rank, stepper_sptr, nu_x, nu_y);
        }
    } else {
        if (rank == 0) {
            std::cout << "===========================" << std::endl;
            std::cout << "   No error element requested" << std::endl;
            std::cout << std::endl;
        }
    }
}

//-----------------------------------------------------------------------------------


Collective_operators get_solver(Cxx_offdiag_options const& opts, int rank)
{

    //Dummy_collective_operator_sptr bpm_measure_sptr;
    //if (opts.bpms){
        //bpm_measure_sptr=Dummy_collective_operator_sptr(new Dummy_collective_operator("bmp_measure"));
    //}

    Collective_operators sc_operators;

    // you can't have multiple space charge solvers enabled
    if ((opts.space_charge_3dh +
         opts.space_charge_2dh +
         opts.space_charge_rec) > 1) {
        throw std::runtime_error("you can't have more than one space charge solver enabled");
    }

    if(opts.solver == "3doh") {
        Space_charge_3d_open_hockney_sptr spc3dh_sptr;
        std::vector<int> grid_shape(3);
        grid_shape[0] = opts.gridx;
        grid_shape[1] = opts.gridy;
        grid_shape[2] = opts.gridz;
        bool longitudinal_kicks = false;
        bool periodic_z = false;
        // We're not using z periodic right now, but need to check this since now most solver
        // stuff is done in the lab frame
        //double z_period=gamma*lattice_length/opts.harmon;
        double z_period = 0.0;
        bool grid_entire_period = false;
        double nsigma=8.;
        Commxx_divider_sptr commxx_divider_sptr= Commxx_divider_sptr(new Commxx_divider(opts.spc_comm_size, false));

        spc3dh_sptr = Space_charge_3d_open_hockney_sptr(new Space_charge_3d_open_hockney(commxx_divider_sptr, grid_shape,longitudinal_kicks,periodic_z,
                    z_period,grid_entire_period,nsigma));
        spc3dh_sptr->set_green_fn_type(Space_charge_3d_open_hockney::linear);
        sc_operators.push_back(spc3dh_sptr);

        // if magic compensation is selected, add an additional compensating space charge kick/
        //   if the compensation degree is larger than 1, you'll have to manually go through the
        //   operators and deactivate the unneeded compensation operators
        if (opts.magiccomp > 0) {
            Space_charge_3d_open_hockney_sptr comp3dh_sptr =
                Space_charge_3d_open_hockney_sptr(
                    new Space_charge_3d_open_hockney(commxx_divider_sptr, grid_shape,longitudinal_kicks,periodic_z,
                        z_period,grid_entire_period,nsigma,-double(opts.magiccomp)));
            sc_operators.push_back(comp3dh_sptr);
        }

        if (rank==0) {
            std::cout<<std::endl;
            std::cout<< "3D HOCKNEY SOLVER WITH OPEN BOUNDARIES"<<std::endl;
            std::cout<<"grid for spc 3dh =["<< grid_shape[0]<<", "
                     << grid_shape[1]<<", "
                     << grid_shape[2]<<"]"
                     <<std::endl;
            std::cout<<"longitudinal kicks ="<<longitudinal_kicks<< std::endl;
            std::cout<<"longitudinal periodicity="<<periodic_z<<std::endl;
            std::cout<<"nsigma="<<nsigma<<std::endl;
            if (opts.magiccomp > 0) {
                std::cout << "magic compensation active at factor: " << opts.magiccomp << std::endl;
            }
            std::cout<<"___________________________________________________________"<<std::endl;
        }
    }  // opts.space_charge_3dh
    else if (opts.solver == "2doh") {
        Space_charge_2d_open_hockney_sptr spc2dh_sptr;
        std::vector<int> grid_shape(3);
        grid_shape[0] = opts.gridx;
        grid_shape[1] = opts.gridy;
        grid_shape[2] = opts.gridz;
        Commxx_divider_sptr commxx_divider_sptr= Commxx_divider_sptr(new Commxx_divider(opts.spc_comm_size, false));
        spc2dh_sptr = Space_charge_2d_open_hockney_sptr(
            new Space_charge_2d_open_hockney(commxx_divider_sptr, grid_shape));
        sc_operators.push_back(spc2dh_sptr);
        if (rank==0) {
            std::cout<<std::endl;
            std::cout<< "2D HOCKNEY SOLVER WITH OPEN BOUNDARIES"<<std::endl;
            std::cout<<"grid for spc 2dh =["<< grid_shape[0]<<", "
                     << grid_shape[1]<<", "
                     << grid_shape[2]<<"]"
                     <<std::endl;
            std::cout<<"___________________________________________________________"<<std::endl;
        }
    } // 2doh
    else if (opts.solver == "dummy") {
        Dummy_collective_operator_sptr dummy_solver_sptr(new Dummy_collective_operator("dummy"));
        if (rank == 0) {
            std::cout << std::endl;
            std::cout << "DUMMY SOLVER NO SPACE CHARGE" << std::endl;
            std::cout<<"___________________________________________________________"<<std::endl;
        }
        sc_operators.push_back(dummy_solver_sptr);
    }
    else {
        throw std::runtime_error("Invalid space charge solver selected: " + opts.solver);
    }
    return sc_operators;
}

//-----------------------------------------------------------------------------------

//-----------------------------------------------------------------------------------

// activate magiccomp
void
activate_magiccomp(Cxx_offdiag_options const& opts, Stepper_sptr stepper_sptr, int rank)
{
    if (opts.magiccomp == 1) {
        // naive compensation at the point of the lens.  This is automatic
        if (rank == 0) {
            std::cout << "Complete compensation activated at all SC elements" << std::endl;
        }
    } else if (opts.magiccomp == 2) {
        if (rank == 0) {
            std::cout << "activating magiccomp 2" << std::endl;
        }
        // naive compensation at factor of two activated.  double every other lens
        // go through the steps, find the space charge operators, every 
        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // if (rank == 0) {
            //     std::cout << "trying step, length:  " << (*sit)->get_length() << ", num operators: " << (*sit)->get_operators().size() << std::endl;
            // }
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        // if (rank == 0) {
                        //     std::cout << "regular SC operator" << std::endl;
                        // }
                        ++scplus;
                    } else if (ks < 0.0) {
                        // if (rank == 0) {
                        //     std::cout << "compensating SC operator" << std::endl;
                        // }
                        ++scminus;
                    }
                }
            }
            // if (rank == 0) {
            //     std::cout << "scplus: " << scplus << ", scminus: " << scminus << std::endl;
            // }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;
            // Now make kick_scale -2 for parity 0 compensating SC operators and remove all parity 1
            // compensating SC operators
            if (sc_parity%2 == 1) {
                // remove the compensating operator from the operator list
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            // if (rank == 0) {
                            //     std::cout << "removing compensating operator" << std::endl;
                            // }
                            (*sit)->get_operators().erase(oit);
                            break;  // finished looping through operators
                        }
                    }
                }
            } else {
                // set the kick strength to -2 for the compensating operator
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            // if (rank == 0) {
                            //     std::cout << "setting strength to -2" << std::endl;
                            // }
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-2.0);
                            break;  // finished looping through operators
                        }
                    }
                }
            }
        }
    }  // end of magiccomp == 2
    else if (opts.magiccomp == 16) {
        if (rank == 0) {
            std::cout << "activating magiccomp 16 (1/6) with compensation of " << opts.sccomp << std::endl;
        }
        // do compensation at strength -1 at 1/6 kicks.  This is meant
        // to compensate the local space charge kick at one step.

        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // if (rank == 0) {
            //     std::cout << "trying step, length:  " << (*sit)->get_length() << ", num operators: " << (*sit)->get_operators().size() << std::endl;
            // }
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "ignoring 0 length step" <<  std::endl;
                }
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "regular SC operator" << std::endl;
                        }
                        ++scplus;
                    } else if (ks < 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "compensating SC operator" << std::endl;
                        }
                        ++scminus;
                    }
                }
            }
            if (DEBUG && (rank == 0)) {
                std::cout << "scplus: " << scplus << ", scminus: " << scminus << std::endl;
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;

            if (DEBUG && (rank == 0)) {
                std::cout << "sc_parity: " << sc_parity << std::endl;;
            }
            // find the compensating kick.  If this is the 1/6
            // step, set the strength to -1, else remove the operator
            // from the operator list

            if (DEBUG && (rank == 0)) {
                std::cout << "loop through operators" << std::endl;
            }
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "operator: " << (*oit)->get_name() << std::endl;
                }
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (DEBUG && (rank == 0)) {
                        std::cout << "space charge 3d open hockney kick strength: " << ks << std::endl;
                    }
                    if (ks < 0.0) {
                        // this is the compensating operator
                        if (sc_parity%6 == 0) {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "setting compensating operator strength to " << -opts.sccomp << std::endl;
                            }
                            // this is the one I keep for compensating, but
                            // its strength is wrong.  Correct it.
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-opts.sccomp);
                            break; // finished with this set of operators, go on to the next step
                        } else {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "removing compensating operator" << std::endl;
                            }
                            (*sit)->get_operators().erase(oit);
                            break; // finished with this set of operators, go on to the next step
                        }
                    }
                }
            }
        }
    } // end of maggiccomp == 16    
    // magiccomp == 13 does compensation at 2/6 of the space charge kick
    // locations
    else if (opts.magiccomp == 13) {
        if (rank == 0) {
            std::cout << "activating magiccomp 13 (2/6) with compensation of " << opts.sccomp << std::endl;
        }
        // do compensation at strength sccomp at 2/6 kicks.

        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // if (rank == 0) {
            //     std::cout << "trying step, length:  " << (*sit)->get_length() << ", num operators: " << (*sit)->get_operators().size() << std::endl;
            // }
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "ignoring 0 length step" <<  std::endl;
                }
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "regular SC operator" << std::endl;
                        }
                        ++scplus;
                    } else if (ks < 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "compensating SC operator" << std::endl;
                        }
                        ++scminus;
                    }
                }
            }
            if (DEBUG && (rank == 0)) {
                std::cout << "scplus: " << scplus << ", scminus: " << scminus << std::endl;
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;

            if (DEBUG && (rank == 0)) {
                std::cout << "sc_parity: " << sc_parity << std::endl;;
            }
            // find the compensating kick.  If this is the 1/6
            // step, set the strength to -1, else remove the operator
            // from the operator list

            if (DEBUG && (rank == 0)) {
                std::cout << "loop through operators" << std::endl;
            }
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "operator: " << (*oit)->get_name() << std::endl;
                }
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (DEBUG && (rank == 0)) {
                        std::cout << "space charge 3d open hockney kick strength: " << ks << std::endl;
                    }
                    if (ks < 0.0) {
                        // this is the compensating operator
                        if ((sc_parity%6 == 0) ||
                            (sc_parity%6 == 3)) {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "setting compensating operator strength to " << -opts.sccomp << std::endl;
                            }
                            // this is the one I keep for compensating, but
                            // its strength is wrong.  Correct it.
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-opts.sccomp);
                            break; // finished with this set of operators, go on to the next step
                        } else {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "removing compensating operator" << std::endl;
                            }
                            (*sit)->get_operators().erase(oit);
                            break; // finished with this set of operators, go on to the next step
                        }
                    }
                }
            }
        }
    } // end of maggiccomp == 13
    // magiccomp == 112 does compensation at 1/12 of the space charge kick
    // locations
    else if (opts.magiccomp == 112) {
        if (rank == 0) {
            std::cout << "activating magiccomp 112 (1/12) with compensation of " << opts.sccomp << std::endl;
        }
        // do compensation at strength sccomp at 2/6 kicks.

        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // if (rank == 0) {
            //     std::cout << "trying step, length:  " << (*sit)->get_length() << ", num operators: " << (*sit)->get_operators().size() << std::endl;
            // }
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "ignoring 0 length step" <<  std::endl;
                }
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "regular SC operator" << std::endl;
                        }
                        ++scplus;
                    } else if (ks < 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "compensating SC operator" << std::endl;
                        }
                        ++scminus;
                    }
                }
            }
            if (DEBUG && (rank == 0)) {
                std::cout << "scplus: " << scplus << ", scminus: " << scminus << std::endl;
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;

            if (DEBUG && (rank == 0)) {
                std::cout << "sc_parity: " << sc_parity << std::endl;;
            }
            // find the compensating kick.  If this is the 1/6
            // step, set the strength to -1, else remove the operator
            // from the operator list

            if (DEBUG && (rank == 0)) {
                std::cout << "loop through operators" << std::endl;
            }
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "operator: " << (*oit)->get_name() << std::endl;
                }
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (DEBUG && (rank == 0)) {
                        std::cout << "space charge 3d open hockney kick strength: " << ks << std::endl;
                    }
                    if (ks < 0.0) {
                        // this is the compensating operator
                        if (sc_parity%12 == 0) {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "setting compensating operator strength to " << -opts.sccomp << std::endl;
                            }
                            // this is the one I keep for compensating, but
                            // its strength is wrong.  Correct it.
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-opts.sccomp);
                            break; // finished with this set of operators, go on to the next step
                        } else {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "removing compensating operator" << std::endl;
                            }
                            (*sit)->get_operators().erase(oit);
                            break; // finished with this set of operators, go on to the next step
                        }
                    }
                }
            }
        }
    } // end of maggiccomp == 112
    // magiccomp == 118 does compensation at 1/18 of the space charge kick
    // locations for a total of 4 compensation kicks/turn
    else if (opts.magiccomp == 118) {
        if (rank == 0) {
            std::cout << "activating magiccomp 112 (1/12) with compensation of " << opts.sccomp << std::endl;
        }
        // do compensation at strength sccomp at 2/18 kicks.

        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // if (rank == 0) {
            //     std::cout << "trying step, length:  " << (*sit)->get_length() << ", num operators: " << (*sit)->get_operators().size() << std::endl;
            // }
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "ignoring 0 length step" <<  std::endl;
                }
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "regular SC operator" << std::endl;
                        }
                        ++scplus;
                    } else if (ks < 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "compensating SC operator" << std::endl;
                        }
                        ++scminus;
                    }
                }
            }
            if (DEBUG && (rank == 0)) {
                std::cout << "scplus: " << scplus << ", scminus: " << scminus << std::endl;
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;

            if (DEBUG && (rank == 0)) {
                std::cout << "sc_parity: " << sc_parity << std::endl;;
            }
            // find the compensating kick.  If this is the 1/6
            // step, set the strength to -1, else remove the operator
            // from the operator list

            if (DEBUG && (rank == 0)) {
                std::cout << "loop through operators" << std::endl;
            }
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "operator: " << (*oit)->get_name() << std::endl;
                }
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (DEBUG && (rank == 0)) {
                        std::cout << "space charge 3d open hockney kick strength: " << ks << std::endl;
                    }
                    if (ks < 0.0) {
                        // this is the compensating operator
                        if (sc_parity%18 == 0) {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "setting compensating operator strength to " << -opts.sccomp << std::endl;
                            }
                            // this is the one I keep for compensating, but
                            // its strength is wrong.  Correct it.
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-opts.sccomp);
                            break; // finished with this set of operators, go on to the next step
                        } else {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "removing compensating operator" << std::endl;
                            }
                            (*sit)->get_operators().erase(oit);
                            break; // finished with this set of operators, go on to the next step
                        }
                    }
                }
            }
        }
    } // end of magiccomp == 118
    else if (opts.magiccomp == 15) {
        if (rank == 0) {
            std::cout << "activating magiccomp 15" << std::endl;
        }
        // eliminate all the compensating kicks.  This should be identical
        // to magiccomp0

        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // if (rank == 0) {
            //     std::cout << "trying step, length:  " << (*sit)->get_length() << ", num operators: " << (*sit)->get_operators().size() << std::endl;
            // }
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                if (DEBUG && (rank == 0)) {
                    std::cout << "ignoring 0 length step" <<  std::endl;
                }
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "regular SC operator" << std::endl;
                        }
                        ++scplus;
                    } else if (ks < 0.0) {
                        if (DEBUG && (rank == 0)) {
                            std::cout << "compensating SC operator" << std::endl;
                        }
                        ++scminus;
                    }
                }
            }
            if (DEBUG && (rank == 0)) {
                std::cout << "scplus: " << scplus << ", scminus: " << scminus << std::endl;
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;
            if (DEBUG && (rank == 0)) {
                std::cout << "sc_parity: " << sc_parity << ": ";
            }

            // eliminate the compensating kick
            if (true) {
                // remove the compensating operator from the operator list
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            if (DEBUG && (rank == 0)) {
                                std::cout << "removing compensating operator" << std::endl;
                            }
                            (*sit)->get_operators().erase(oit);
                            break;  // finished looping through operators
                        }
                    }
                }
            }
        }
    }   // end of maggiccomp == 15
    else if (opts.magiccomp == 3) {
        if (rank == 0) {
            std::cout << "activating magiccomp: " << opts.magiccomp << std::endl;
        }
        // do the same thing as magiccomp == 2 except opposite parity
        // naive compensation at factor of two activated.  double every other lens
        // go through the steps, find the space charge operators, every 
        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        ++scplus;
                    } else if (ks < 0.0) {
                        ++scminus;
                    }
                }
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;
            // Now make kick_scale -2 for parity 1 compensating SC operators and remove all parity 0
            // compensating SC operators
            // if (rank == 0) {
            //     std::cout << "sc_parity: " << sc_parity << ", sc_parity%2: " << sc_parity%2 << std::endl;
            // }
            if (sc_parity%2 == 0) {
                // if (rank == 0) {
                //     std::cout << "removing operator" << std::endl;
                // }
              // remove the compensating operator from the operator list
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            (*sit)->get_operators().erase(oit);
                            break;  // finished looping through operators
                        }
                    }
                }
            } else {
                // set the kick strength to -2 for the compensating operator
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-2.0);
                            break;  // finished looping through operators
                        }
                    }
                }
            }
        }
    } else if (opts.magiccomp == 4) {
        // Use lens2 to compensate for space charge in first half cell.  Scale the kick to account for phase advance.
        // The kick is scaled to compensate for large action 0.
        const double phi_1_2 = 1.188194066316849 -  0.7700655111169677;
        const double alpha2 = -1.843930877445086;
        const double F = (-alpha2 * std::tan(phi_1_2) + 1.0 );
        if (rank == 0) {
            std::cout << "magiccomp 4 scaling kick strength by F = cos(phi)*(-alpha_2 * sin(phi) + cos(phi))" << std::endl;
            std::cout << "F: " << F << std::endl;
        }
        
        int sc_parity = 0;
        int stepcnt = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit, ++stepcnt) {
            // if (rank == 0) {
            //     std::cout << "\nnew step " << stepcnt << ", length:  " << (*sit)->get_length() << ", num operators: " << (*sit)->get_operators().size() << std::endl;
            // }
            // print_sc_scales(stepper_sptr);
            // (*sit)->print(stepcnt);

            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                //(*oit)->print();
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        // if (rank == 0) {
                        //     std::cout << "regular SC operator" << std::endl;
                        // }
                        ++scplus;
                    } else if (ks < 0.0) {
                        // std::cout << "compensating SC operator" << std::endl;
                        ++scminus;
                    }
                }
            }
            if (rank == 0) {
                // std::cout << "scplus: " << scplus << ", scminus: " << scminus << std::endl;
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;

            // std::cout << "sc_parity: " << sc_parity << ", sc_parity%2: " << sc_parity%2 << std::endl;
            
            // parity 1 steps contain lens1. Deactivate the conpensating operator.
            if (sc_parity%2 == 1) {
                // if (rank == 0) {
                //     std::cout << "this is parity 1, removing compensating operator" << std::endl;
                // }
                // remove the compensating operator from the operator list
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    // std::cout << "operator: " << (*oit)->get_name() << (*oit)->get_type() << std::endl;
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            // std::cout << "removing compensating operator" << std::endl;
                            (*sit)->get_operators().erase(oit);
                            // std::cout << "setting compensating operator strength to 0" << std::endl;
                            // boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(0.0);
                            break;  // finished looping through operators
                        }
                    }
                }
            } else {
                // if (rank == 0) {
                //     std::cout << "this is parity 0, setting strength of compensating operator" << std::endl;
                // }
                // scale the kick strength to by factor -F for the compensating operator
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    // std::cout << "operator: " << (*oit)->get_name() << (*oit)->get_type() << std::endl;
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            // std::cout << "setting strength of operator" << std::endl;
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-F-1.0);
                            break;  // finished looping through operators
                        }
                    }
                }
            }
        }
    } else if (opts.magiccomp == 5) {
        // Use lens2 to compensate for space charge in first half cell.  Scale the kick to account for phase advance.
        // The kick is scaled to compensate for action 0.
        const double phi_1_2 = 1.188194066316849 -  0.7700655111169677;
        const double alpha2 = -1.843930877445086;
        const double F =  (-alpha2 * std::tan(phi_1_2) + 1.0 ) * std::cos(phi_1_2)*std::cos(phi_1_2);
        if (rank == 0) {
            std::cout << "magiccomp 5 scaling kick strength by F = cos(phi)*(-alpha_2 * sin(phi) + cos(phi)) * cos(phi)**2" << std::endl;
            std::cout << "F: " << F << std::endl;
        }
        
        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        ++scplus;
                    } else if (ks < 0.0) {
                        ++scminus;
                    }
                }
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;
            // parity 1 steps contain lens1. Deactivate the conpensating operator.
            if (sc_parity%2 == 1) {
                // remove the compensating operator from the operator list
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            (*sit)->get_operators().erase(oit);
                            break;  // finished looping through operators
                        }
                    }
                }
            } else {
                // scale the kick strength to by factor F for the compensating operator
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-F-1.0);
                            break;  // finished looping through operators
                        }
                    }
                }
            }
        }
    } else if (opts.magiccomp == 6) {
        // compensate perfectly for half of the lenses by removing
        // half the compensation
        if (rank == 0) {
            std::cout << "magiccomp 6  removing half the compensation lenses";
        }
        
        int sc_parity = 0;
        for (Steps::iterator sit=stepper_sptr->get_steps().begin(); sit!=stepper_sptr->get_steps().end(); ++sit) {
            // 0 length step can't have anything useful in it
            if ((*sit)->get_length() == 0.0) {
                continue;
            }
            // this step better have 2 space charge 3d operators.  One with positive kick_scale and one with negative kick_scale
            int scplus=0;
            int scminus=0;
            for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                if ((*oit)->get_name() == "space charge 3D open hockney") {
                    int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                    if (ks > 0.0) {
                        ++scplus;
                    } else if (ks < 0.0) {
                        ++scminus;
                    }
                }
            }
            if ((scplus != 1) || (scminus != 1)) {
                if (rank == 0) {
                    std::cout << "skipping step without the right number of SC3DOH operators" << std::endl;
                }
                continue; // go to next step
            }
            ++sc_parity;
            // parity 1 steps contain lens1. Deactivate the conpensating operator.
            if (sc_parity%2 == 1) {
                // remove the compensating operator from the operator list
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            (*sit)->get_operators().erase(oit);
                            break;  // finished looping through operators
                        }
                    }
                }
            } else {
                // scale the kick strength to by factor 1 for total compensation
                for (Operators::iterator oit=(*sit)->get_operators().begin(); oit!=(*sit)->get_operators().end(); ++oit) {
                    if ((*oit)->get_name() == "space charge 3D open hockney") {
                        int ks = boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->get_kick_scale();
                        if (ks < 0.0) {
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*oit)->set_kick_scale(-1.0);
                            break;  // finished looping through operators
                        }
                    }
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------------
void
run(Cxx_offdiag_options const& opts)
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // std::cout<<"run start on rank: "<<rank<<std::endl;

    Lattice_sptr lattice_sptr(setup_lattice(opts, rank));


    MPI_Barrier(MPI_COMM_WORLD);
    double tini = MPI_Wtime();


    Reference_particle reference_particle=lattice_sptr->get_reference_particle();
    double lattice_length=lattice_sptr->get_length();
    double beta = reference_particle.get_beta();
    double gamma = reference_particle.get_gamma();
    double energy = reference_particle.get_total_energy();


    if (rank==0) {
        std::cout<<"    ***     before stepper   ***     "<<std::endl;
        std::cout<<" beta="<<beta<<std::endl;
        std::cout<<" gamma="<<gamma<<std::endl;
        std::cout<<" lattice_length="<<lattice_length<<std::endl;
        std::cout<<" energy="<<energy<<std::endl;
        std::cout<<std::endl;
        std::cout<<"    ***    ***********    ***     "<<std::endl;
    }

    Stepper_sptr stepper_sptr;

    if (opts.stepper == "independent") {
        stepper_sptr = Stepper_sptr(new Independent_stepper(lattice_sptr, opts.map_order, opts.num_steps));
        if (rank==0) {
            std::cout<<std::endl;
            std::cout<<"no collective effects, no bpm, independent stepper with " << opts.num_steps << " steps" <<std::endl;
            std::cout<<"___________________________________________________________"<<std::endl;
        }
    } else if (opts.stepper == "elements")  {
        stepper_sptr = Stepper_sptr(new Independent_stepper_elements(lattice_sptr, opts.map_order, opts.num_steps));
        if (rank==0) {
            std::cout<<std::endl;
            std::cout<<"no collective effects, no bpm, independent element stepper with " << opts.num_steps << " steps/element" <<std::endl;
            std::cout<<"___________________________________________________________"<<std::endl;
        }
    } else if (opts.stepper == "splitoperator") {
        if (rank == 0) {
            std::cout<<"split operator stepper with " << opts.num_steps << " steps/turn" <<std::endl;
            std::cout<<"___________________________________________________________"<<std::endl;
        }
        //if (opts.space_charge_3dh) {
            //operators.push_back(spc3dh_sptr);
        //} else {
            //operators.push_back(Dummy_collective_operator_sptr(new Dummy_collective_operator("dummy_spc")));
        //}
        stepper_sptr = Stepper_sptr(new Split_operator_stepper(lattice_sptr, opts.map_order, get_solver(opts, rank), opts.num_steps));
    } else if (opts.stepper == "splitoperatorelements") {
        if (rank == 0) {
            std::cout<<"split operator stepper elements with " << opts.num_steps << " steps/element" <<std::endl;
            std::cout<<"___________________________________________________________"<<std::endl;
        }
        stepper_sptr = Stepper_sptr(new Split_operator_stepper_elements(lattice_sptr, opts.map_order, get_solver(opts, rank), opts.num_steps));
    }
#if 0
    else if (opts.stepper == "choice") {
        // if ( (opts.space_charge_2dh) || (opts.space_charge_3dh) || (opts.space_charge_rec) || (opts.bpms) ) {
        List_choice_map list_choice_map;

        Collective_operators operators_quad;

        if (opts.space_charge_3dh) operators_quad.push_back(spc3dh_sptr);
        Kicks  kicks_quad(operators_quad,opts.steps_per_quad);

        int steps_per_bpm=1;
        Collective_operators operators_bpm;
        if (opts.bpms) operators_bpm.push_back(bpm_measure_sptr);
        Kicks  kicks_bpm(operators_bpm, steps_per_bpm);

        for (Lattice_elements::const_iterator latt_it =
                 lattice_sptr->get_elements().begin();
             latt_it !=lattice_sptr->get_elements().end(); ++latt_it){
            if (opts.space_charge_3dh){
                if (((*latt_it)->get_name()=="f") || ((*latt_it)->get_name()=="d")) {
                    // std::cout<<" element_name: "<<(*latt_it)->get_name()<<std::endl;
                    list_choice_map[(*latt_it)->get_name()]=kicks_quad;
                }
            }
            if (opts.bpms) {
                if((*latt_it)->get_name()=="bpm") {
                    list_choice_map[(*latt_it)->get_name()]=kicks_bpm;
                    // if (rank==0) std::cout<<" BPMS are: "<<(*latt_it)->get_name()<<std::endl;
                }
            }


        }

        Collective_operators operators_else;
        if (opts.space_charge_3dh) operators_else.push_back(spc3dh_sptr);
        Kicks  kicks_else(operators_else,opts.num_steps_else);
        list_choice_map["else"]=kicks_else;

        stepper_sptr=Stepper_sptr(new Split_operator_stepper_choice(lattice_sptr, opts.map_order, list_choice_map, true));
        if (rank==0){
            std::cout<<std::endl;
            std::cout<<"Split_operator_stepper_choice created"<<std::endl;
            if (opts.bpms) std::cout<<"measurements at BPMS elements"<<std::endl;
            std::cout<<"steps_per_quad="<<opts.steps_per_quad<<std::endl;
            std::cout<<"num_steps_else="<<opts.num_steps_else<<std::endl;
            std::cout<<"___________________________________________________________"<<std::endl;
        }
    }
#endif
    else {
        throw std::runtime_error("unknown stepper choice: " + opts.stepper);
    }

    std::pair<double,double> orig_tunes = stepper_sptr->get_lattice_simulator().get_both_tunes(true);
    if (rank == 0) {
        std::cout << "Original lattice tunes: hor: " << orig_tunes.first << ", ver: " << orig_tunes.second << std::endl;
    }

    double new_nux = orig_tunes.first;
    double new_nuy = orig_tunes.second;
    bool need_tune_adjust = false;
    if (opts.set_nux != 0.0) {
        new_nux = opts.set_nux;
        need_tune_adjust = true;
    }
    if (opts.set_nuy != 0.0) {
        new_nuy = opts.set_nuy;
        need_tune_adjust = true;
    }
    if (need_tune_adjust) {
        set_tunes(rank, stepper_sptr, new_nux, new_nuy);
        if (rank == 0) {
            std::cout << "Setting tunes to: " << new_nux << ", " << new_nuy << std::endl;
        }
    }

    // activate magic compensation                                             
    activate_magiccomp(opts, stepper_sptr, rank);

    // activate error element
    activate_error_element(opts, rank, stepper_sptr, new_nux, new_nuy);

    stepper_sptr->get_lattice_simulator().register_closed_orbit();
    stepper_sptr->get_lattice_simulator().set_rf_bucket_length();
    if (rank == 0) {
        stepper_sptr->get_lattice_simulator().print_lattice_functions();
    }

    if (rank==0) {
        std::cout<<"stepper:lattice_simulator: map order="<< stepper_sptr->get_lattice_simulator().get_map_order() <<std::endl;
        std::cout<<"stepper: number of steps="<<stepper_sptr->get_steps().size()<<std::endl;
        std::cout<<"stepper rf frequency="<<stepper_sptr->get_lattice_simulator().get_rf_frequency()<<std::endl;
        //stepper_sptr->print();
    }

    if (rank == 0) {
        print_lattice(stepper_sptr->get_lattice_simulator().get_lattice_sptr(), "adjusted_lattice.txt");
        write_lsx_lattice(stepper_sptr->get_lattice_simulator().get_lattice_sptr(), "adjusted_lattice.lsx");
    }

    stepper_sptr->print_cs_step_betas();
    reference_particle=stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_reference_particle();
    lattice_length=stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_length();
    beta = reference_particle.get_beta();
    gamma = reference_particle.get_gamma();
    energy = reference_particle.get_total_energy();
    //  double bunch_sp=stepper_sptr->get_lattice_simulator().get_bucket_length();
    MArray1d clo=stepper_sptr->get_lattice_simulator().get_closed_orbit();

    if (rank==0) {
        std::cout<<"    ********    stepper  lattice  ************     "<<std::endl;
        std::cout<<std::endl;
        multi_array_print(clo,"closed_orbit");
        //        std::cout<<" bunch spacing="<<bunch_sp<<std::endl;
        std::cout<<" rf bucket length="<<stepper_sptr->get_lattice_simulator().get_bucket_length()<<" ="
                 <<stepper_sptr->get_lattice_simulator().get_bucket_length()*1e9/beta/pconstants::c<<" ns"<<std::endl;
        std::cout<<" rf frequency="<<stepper_sptr->get_lattice_simulator().get_rf_frequency()/1.e6<<" MHz"<<std::endl;
        std::cout<<" beta="<<beta<<std::endl;
        std::cout<<" gamma="<<gamma<<std::endl;
        std::cout<<" lattice_length="<<lattice_length<<std::endl;
        std::cout<<" closed_orbit_length="<<stepper_sptr->get_lattice_simulator().get_closed_orbit_length()<<std::endl;
        std::cout<<" energy="<<energy<<"GeV"<<std::endl;
        std::cout<<" reference momentum="<<reference_particle.get_momentum()<<" GeV/c"<<std::endl;
        std::cout<<std::endl;
        std::cout<<"    ***********************************     "<<std::endl;
        std::cout<<std::endl;
    }

    // print out the lattice functions at the boc points
    if (rank == 0) {
        std::ofstream lfo("boc_lf.txt", std::ofstream::out);
        std::cout << "\nLattice functions at boc and lenses" << std::endl;
        Lattice_sptr lattice_sptr(stepper_sptr->get_lattice_simulator().get_lattice_sptr());
        for (Lattice_elements::const_iterator it = lattice_sptr->get_elements().begin();
            it != lattice_sptr->get_elements().end(); ++it) {

            if ( ((*it)->get_name() == "boc") && ((*it)->get_type() == "marker") ||
                ((*it)->get_name() == "lens1") && ((*it)->get_type() == "elens") ||
                ((*it)->get_name() == "lens2") && ((*it)->get_type() == "elens") ) {
                Lattice_functions lf(stepper_sptr->get_lattice_simulator().get_lattice_functions(**it));
                std::cout << lf.arc_length << ": " << (*it)->get_name() << " : betax: " << lf.beta_x << ", alphax: " << lf.alpha_x << ", betay: " << lf.beta_y << ", alphay: " << lf.alpha_y << std::endl;
                if ((*it)->get_name() == "boc") {
                    lfo << lf.arc_length << ":" << lf.beta_x << ":" << lf.alpha_x << ":" << lf.beta_y << ":" << lf.alpha_y << std::endl;
                }
            }
        }
        lfo.close();
    }

    // these are the lattice functions at the beginning of the line
    Lattice_functions lf(stepper_sptr->get_lattice_simulator().get_lattice_functions(*lattice_sptr->get_elements().back()));

    double emitx = opts.emitx;
    double emity = opts.emity;
    double stdx = std::sqrt(emitx*lf.beta_x);
    double stdy = std::sqrt(emity*lf.beta_y);
    double mean_sigma = std::sqrt(stdx*stdy);
    double aperture_radius = opts.aperture_sigma*mean_sigma;

    if (rank == 0) {
        std::cout << "X emittance: " << emitx << std::endl;
        std::cout << "stdx: " << stdx << std::endl;
        std::cout << "Y emittance: " << emity << std::endl;
        std::cout << "stdy: " << stdy << std::endl;
        std::cout << "tune x: " << lf.psi_x/(2.0*mconstants::pi) << std::endl;
        std::cout << "tune y: " << lf.psi_y/(2.0*mconstants::pi) << std::endl;
    }

    MArray2d correlation_matrix(boost::extents[6][6]);
    if (opts.map_order==1){
        MArray2d one_turn_map=stepper_sptr->get_lattice_simulator().get_linear_one_turn_map();
        if (opts.flatbucket) {
            std::vector <int> rms_index(3);
            rms_index[0] = Bunch::x;
            rms_index[1] = Bunch::y;
            rms_index[2] = Bunch::zp;
            correlation_matrix = get_correlation_matrix(one_turn_map, stdx, stdy, opts.dpoprms, beta, rms_index);
        } else {
            correlation_matrix=get_correlation_matrix(one_turn_map,stdx, stdy, opts.zrms, beta);
        }
        if (opts.beamparams) {
            if (rank == 0) {
                double emitx = opts.transhemit;
                double emity = opts.transvemit;
                double betax = opts.hbeta;
                double betay = opts.vbeta;
                double alphax = opts.halpha;
                double alphay = opts.valpha;
                
                std::cout << "using transverse parameters for bunch:" << std::endl;
                std::cout << "horiz emit: " << emitx << ", vert emit" << emity << std::endl;
                std::cout << "horiz beta: " << betax << ", vert beta: " << betay << std::endl;
                std::cout << "horiz alpha: " << alphax << ", vert alpha: " << alphay << std::endl;
                for (int i=0; i<4; ++i) {
                    for (int j=0; j<4; ++j) {
                        correlation_matrix[i][j] = 0.0;
                    }
                }
                for (int i=4; i<6; ++i) {
                    for (int j=0; j<4; ++j) {
                        correlation_matrix[i][j] = 0.0;
                        correlation_matrix[j][i] = 0.0;
                    }
                }
                correlation_matrix[0][0] = emitx*betax;
                correlation_matrix[0][1] = -emitx*alphax;
                correlation_matrix[1][0] = -emitx*alphax;
                correlation_matrix[1][1] = emitx * (1 + alphax*alphax)/betax;
                correlation_matrix[2][2] = emity*betay;
                correlation_matrix[2][3] = -emity*alphay;
                correlation_matrix[3][2] = -emity*alphay;
                correlation_matrix[3][3] = emity * (1 + alphay*alphay)/betay;
            }
        }
        if (rank==0){
            std::cout<<" correlation matrix="<<multi_array_to_string(correlation_matrix)<<std::endl;
        }
        int map_size=one_turn_map.size();
        Eigen::MatrixXd eigen_map(map_size,map_size);
        for (int i=0;i<map_size;++i){
            for (int j=0;j<map_size;++j){
                eigen_map(i,j)=one_turn_map[i][j];
            }
        }
        if (rank==0){
            std::cout<<" one turn map="<<multi_array_to_string(one_turn_map)<<std::endl;
            std::cout<<"eigenvalues of the one_turn_map:"
                     <<"\n"<<"**********************"<<std::endl;
            for (int i=0;i<map_size;++i){
                std::complex<double> vv=eigen_map.eigenvalues()(i);
                std::cout<<"eigenvalue["<<i<<"]="<<vv
                         <<",  absolute value="<<abs(vv)
                         <<",  fractional tune="<<fabs(log(vv).imag()/(2.*mconstants::pi))
                         <<",  "<<1.-fabs(log(vv).imag()/(2.*mconstants::pi))<<std::endl;
            }
            std::cout<<"**********************"<<std::endl;
        }
    }
    double chef_frac_tunex=stepper_sptr->get_lattice_simulator().get_horizontal_tune();
    double chef_frac_tuney=stepper_sptr->get_lattice_simulator().get_vertical_tune();
    double chef_eigen_tunex=stepper_sptr->get_lattice_simulator().get_horizontal_tune(1);
    double chef_eigen_tuney=stepper_sptr->get_lattice_simulator().get_vertical_tune(1);
    double horizontal_chromaticity=stepper_sptr->get_lattice_simulator().get_horizontal_chromaticity();
    double vertical_chromaticity=stepper_sptr->get_lattice_simulator().get_vertical_chromaticity();
    double momentum_compaction=stepper_sptr->get_lattice_simulator().get_momentum_compaction();
    double slip_factor=stepper_sptr->get_lattice_simulator().get_slip_factor();


    if (rank==0){
        std::cout<< "chef FracTune x: "<< chef_frac_tunex<< ", EigenTune x: "<< chef_eigen_tunex<<std::endl;
        std::cout<< "chef FracTune y: "<< chef_frac_tuney<< ", EigenTune y: "<< chef_eigen_tuney<<std::endl;
        std::cout<< "horizontal chromaticity: "<< horizontal_chromaticity<<std::endl;
        std::cout<< "vertical   chromaticity: "<< vertical_chromaticity<<std::endl;
        std::cout<< "momentum compaction: "<< momentum_compaction<<std::endl;
        std::cout<< "slip factor: "<< slip_factor<<std::endl;
    }

    Bunches bunches;
    Commxx_sptr parent_comm_sptr(new Commxx);
    Commxxs comms(generate_subcomms(parent_comm_sptr, opts.num_bunches));
    const double mpweight = opts.real_particles/opts.macroparticles;
    for (int i=0;i<opts.num_bunches;++i) {
        Commxx_sptr commx=comms[i];
        // if zparticle is true, generate a bunch with a single particle at
        // 0.  The rest will be filled in later.
        Bunch_sptr bunch_sptr;
        if (opts.zparticle) {
            bunch_sptr.reset(new Bunch(stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle(),
                                                  1, mpweight, commx));
            MArray2d_ref local_particles(bunch_sptr->get_local_particles());
            if (commx->get_rank() == commx->get_size()-1) {
                for (int i=0; i<6; ++i) {
                    local_particles[0][i] = 0.0;
                }
            }
        } else { // ! opts.zparticle
            bunch_sptr.reset(new Bunch(stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle(),
                                                  opts.macroparticles, opts.real_particles, commx));
        }
        bunch_sptr->set_bucket_index(i);
        bunch_sptr->set_z_period_length(lattice_sptr->get_length()/opts.harmon) ;
        if (commx->has_this_rank()){
            if (opts.load_bunch){
                bunch_sptr->read_file("el_input_particles.h5");
                if (commx->get_rank()==0) std::cout<<" bunch "<<i<<" loaded from:  "<<"el_input_particles.h5"<<std::endl;
            }
            else if (!opts.flatbucket) { // gaussian bucket
                // replace this with truncated normal form
#if 1
                // populate gaussian distribution wwith correlation matrix
                if (rank==0) {
                    std::cout<<" input (xrms, yrms, zrms)[m] = ("<< stdx<<", "<< stdy<<", "<< opts.dpoprms << ")" << std::endl;
                }
                Random_distribution dist(opts.seed,*commx);
                MArray1d input_means(boost::extents[6]);
                for(int imean=0;imean<6;++imean){
                    input_means[imean]=0.;
                }
                if (opts.map_order==1){
                    if (!opts.zparticle) {
                        populate_6d(dist, *bunch_sptr, input_means, correlation_matrix);
                    } else {
                        Bunch_sptr newbunch_sptr(new Bunch(stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle(), (opts.macroparticles-1), mpweight*(opts.macroparticles-1), commx));
                        populate_6d(dist, *newbunch_sptr, input_means, correlation_matrix);
                        bunch_sptr->inject(*newbunch_sptr);
                    }
                }
                else{
                    throw  std::runtime_error("bunch populate for high order maps nor implemented yet");
                }
#else
                // generate with truncated gaussian from normal forms
                Random_distribution dist(opts.seed,*commx);
                MArray1d limits(boost::extents[3]);
                limits[0] = 10.0;
                limits[1] = 10.0;
                limits[2] = stepper_sptr->get_lattice_simulator().get_bucket_length()/opts.zrms;
                MArray1d input_means(boost::extents[6]);
                for(int imean=0;imean<6;++imean){
                    input_means[imean]=0.;
                }
                std::vector<double >  actions(3);
                actions= stepper_sptr->get_lattice_simulator().get_stationary_actions(stdx, stdy, opts.zrms/beta);
                 populate_6d_stationary_gaussian_truncated(dist, *bunch_sptr,  actions, stepper_sptr->get_lattice_simulator(),limits);
#endif
            } // gaussian bucket (end of !flatbucket)
            else if (opts.flatbucket) {  // flatbucket

                if (rank==0) {
                    std::cout<<" longitudinally flat bucket (xrms, yrms, dpoprms)[m] = ("<< opts.xrms<<", "<< opts.yrms<<", " << opts.dpoprms << ")"<<std::endl;
                }
                Random_distribution dist(opts.seed,*commx);
                MArray1d input_means(boost::extents[6]);
                for(int imean=0;imean<6;++imean){
                    input_means[imean]=0.;
                }
                // generate the transverse and dp/p coordinates
                populate_6d(dist, *bunch_sptr, input_means, correlation_matrix);
                // populate the longitudinal coordinate uniformly
                double bucket_length = stepper_sptr->get_lattice_simulator().get_bucket_length();
                populate_longitudinal_uniform(dist, *bunch_sptr, bucket_length);
                // turn on bucket barrier conditions
                bunch_sptr->set_bucket_barrier_length(bucket_length);
                // turn off the RF cavities for bucket barrier
                for (Lattice_elements::iterator leit=stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_elements().begin();
                     leit != stepper_sptr->get_lattice_simulator().get_lattice_sptr()->get_elements().end(); ++leit) {
                    if ((*leit)->get_type() == "rfcavity") {
                        (*leit)->set_double_attribute("volt", 0.0);
                    }
                }
                stepper_sptr->get_lattice_simulator().update();
            }  // end of flatbucket
            else {
                throw std::runtime_error("unknown bucket shape error");
            }
        } //has_this_rank
        bunches.push_back(bunch_sptr);
    }// i

    Bunch_train_sptr bunch_train_sptr(new Bunch_train(bunches, stepper_sptr->get_lattice_simulator().get_bucket_length()));
    Bunch_train_simulator bunch_train_simulator(bunch_train_sptr);

    std::list<int> step_numbers;
    int step_number=0;

    for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it != stepper_sptr->get_steps().end(); ++it) {
        ++step_number;
        if (opts.bpms) {
            for (Operators::const_iterator o_it = (*it)->get_operators().begin(); o_it
                     != (*it)->get_operators().end(); ++o_it) {
                //  std::cout<<"(*o_it)->get_name()="<< (*o_it)->get_name()<<std::endl;
                if ((*o_it)->get_name()=="bmp_measure") step_numbers.push_back(step_number);
            }
        }
        else {
            step_numbers.push_back(step_number);
        }
    }
#if 0
    if (rank == 0) {
        std::cout << " number of measurements per turn="<<step_numbers.size()<<std::endl;
        std::cout << "  saving tracks at steps: ";
        for (std::list<int>::const_iterator it=step_numbers.begin(); it!=step_numbers.end(); ++it) {
            if (it != step_numbers.begin()) std::cout << ", ";
            std::cout << (*it);
        }
        std::cout << std::endl;
    }
#endif

    for (int i=0; i<bunch_train_sptr->get_size(); ++i) {
        std::stringstream bunch_label;
        bunch_label<<i;

#define DIAGNOSTICS 1

#if DIAGNOSTICS
        // loss diagnostics
        if (bunch_train_sptr->get_bunches()[0]->get_comm().has_this_rank()) {
            Diagnostics_loss_sptr diag_loss_sptr(new Diagnostics_loss("apertures_loss.h5", "aperture"));
            diag_loss_sptr->set_bunch_sptr(bunch_train_sptr->get_bunches()[0]);
            stepper_sptr->get_lattice_simulator().get_lattice_sptr()->add_loss_diagnostics(diag_loss_sptr);
        }

        Diagnostics_sptr diagnostics_full2_sptr(new Diagnostics_full2("full2_"+bunch_label.str()+".h5"));
        bunch_train_simulator.add_per_turn(i, diagnostics_full2_sptr);

        Diagnostics_sptr cell_diagnostics_full2_sptr(new Diagnostics_full2("cell_full2_"+bunch_label.str()+".h5"));
        // this is old
        bunch_train_simulator.add_per_forced_diagnostics_step(i, cell_diagnostics_full2_sptr);
#endif

#if 0
        // this is old
        bunch_train_simulator.add_per_step(i,diagnostics_full2_sptr, step_numbers );
#endif

#if DIAGNOSTICS
        if (opts.turn_particles) {
            Diagnostics_sptr diagnostics_particles_sptr(new Diagnostics_particles("turn_particles_"+bunch_label.str()+".h5"));
            bunch_train_simulator.add_per_turn(i,diagnostics_particles_sptr,opts.turn_period);
        }

        if (opts.cell_particles) {
            Diagnostics_sptr cell_diagnostics_particles_sptr(new Diagnostics_particles("cell_particles_"+bunch_label.str()+".h5"));
            bunch_train_simulator.add_per_forced_diagnostics_step(i, cell_diagnostics_particles_sptr);
        }

        if(opts.turn_track){
            int nsaved=std::min(opts.turn_track, opts.macroparticles);
            Diagnostics_sptr diagnostics_bulk_track_sptr(new Diagnostics_bulk_track("bulk_track_"+bunch_label.str()+".h5", nsaved ));
#if 0
            bunch_train_simulator.add_per_step(i,diagnostics_bulk_track_sptr, step_numbers );
#endif
            bunch_train_simulator.add_per_forced_diagnostics_step(i,diagnostics_bulk_track_sptr);
        }
        //             if(opts.phase_space){
        //                 int grid_z=20;
        //                 int grid_zp=20;
        //                 double z_nsigma=4.0;
        //                 double zp_nsigma=4.0;
        //                 Diagnostics_sptr diagnostics_phase_space_density_sptr(new
        //                         Diagnostics_phase_space_density("phase_space_density_"+bunch_label.str()+".h5",
        //                                                         grid_z,grid_zp,z_nsigma,zp_nsigma));
        //                 bunch_train_simulator.add_per_turn(i, diagnostics_phase_space_density_sptr);
        //             }
        //


        if (opts.space_charge_3dh && opts.spc_tuneshift) {
            Diagnostics_space_charge_3d_hockney_sptr sp_diagnostics_sptr
                (new Diagnostics_space_charge_3d_hockney("space_charge_3dh_diagnostics_"+bunch_label.str()+".h5"));
            sp_diagnostics_sptr->set_bunch_sptr(bunch_train_simulator.get_bunch_train().get_bunches()[i]);

            // Add these diagnostics to all the space charge operators
            for (Steps::const_iterator it = stepper_sptr->get_steps().begin(); it != stepper_sptr->get_steps().end(); ++it) {
                //std::cout << "step: length: " << (*it)->get_length() << std::endl;
                for (Operators::const_iterator o_it = (*it)->get_operators().begin(); o_it
                    != (*it)->get_operators().end(); ++o_it) {
                    bool added_one = false;
                    if ((*o_it)->get_name()=="space charge 3D open hockney") {
                        //std::cout << "space charge 3d open hockney kick_scale: " << boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*o_it)->get_kick_scale() << std::endl;
                        //std::cout << "before adding diagnostics " << boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*o_it)->has_diagnostics() << std::endl;
                        //boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*o_it)->add_diagnostics(sp_diagnostics_sptr);
                        if (added_one || !boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*o_it)->has_diagnostics()) {
                            boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*o_it)->add_diagnostics(sp_diagnostics_sptr);
                            added_one = true;
                            break;
                        }
                        //std::cout << "after adding diagnostics " << boost::dynamic_pointer_cast<Space_charge_3d_open_hockney>(*o_it)->has_diagnostics() << std::endl;
                    }
                }
            }
        }
#endif // DIAGNOSTICS

        //adjust means
        if (bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm().has_this_rank()){
            MArray1d bunch_means=Core_diagnostics::calculate_mean(*bunch_train_simulator.get_bunch_train().get_bunches()[i]);
            MArray2d_ref particles( bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_local_particles());
            for (int part=0;part<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_local_num();++part){
                particles[part][0] +=-bunch_means[0]+opts.xoffset;//+clo[0];
                particles[part][1] +=-bunch_means[1];//+clo[1];
                particles[part][2] +=-bunch_means[2]+opts.yoffset;//+clo[2];
                particles[part][3] +=-bunch_means[3];//+clo[3];
                particles[part][4] +=-bunch_means[4]+opts.zoffset;//+beta*clo[4];
                particles[part][5] +=-bunch_means[5];//+clo[5];
            }
            bunch_means=Core_diagnostics::calculate_mean(*bunch_train_simulator.get_bunch_train().get_bunches()[i]);
            MArray1d bunch_stds=Core_diagnostics::calculate_std(*bunch_train_simulator.get_bunch_train().get_bunches()[i],bunch_means);
            MArray2d bunch_mom2=Core_diagnostics::calculate_mom2(*bunch_train_simulator.get_bunch_train().get_bunches()[i],bunch_means);
            if (bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_comm().get_rank()==0) {
                std::cout<<std::endl;
                std::cout<<"bunch # "<<i<<" is perriodic="<<bunch_train_simulator.get_bunch_train().get_bunches()[i]->is_z_periodic()<<std::endl;
                std::cout<<"bunch # "<<i<<" has longitudinal aperture="<<bunch_train_simulator.get_bunch_train().get_bunches()[i]->has_longitudinal_aperture()<<std::endl;
                std::cout<<"bunch # "<<i<<" number of real  particles= "
                         <<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_real_num()<<std::endl;
                std::cout<<"bunch # "<<i<<" number of macroparticles= "
                         <<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_total_num()<<std::endl;
                std::cout<<"bunch # "<<i<<" bucket index= "
                         <<bunch_train_simulator.get_bunch_train().get_bunches()[i]->get_bucket_index()<<std::endl;
                std::cout<<"bunch # "<<i<<" initial offsets (x,xp,y,yp,ct,dpp)=("
                         <<bunch_means[0]<<", "<<bunch_means[1]<<", "<<bunch_means[2]<<", "
                         <<bunch_means[3]<<", "<<bunch_means[4]<<", "<<bunch_means[5]
                         <<") [meters]"<<std::endl;
                std::cout<<"bunch # "<<i<<" initial stds (xrms,yrms,zrms)=("
                         <<bunch_stds[0]<<", "<<bunch_stds[2]<<", "<<beta*bunch_stds[4]<<") [meters]"<<std::endl;
                std::cout<<"bunch # "<<i<<" :"<<std::endl;
                Core_diagnostics::print_bunch_parameters(bunch_mom2, beta);
                std::cout<<"___________________________________________________________"<<std::endl;

            }
        } // has rank

    }// for i

    if (rank == 0) {
        std::cout << "apertures are: ";
        if (!opts.apertures) {
            std::cout << "OFF" << std::endl;;
        } else {
            std::cout << "ON" << std::endl;;
            std::cout << "    aperture radius: " << opts.aperture_sigma << " sigma" << std::endl;
        }
    }
    if (opts.apertures) {
        activate_apertures(opts, stepper_sptr, rank, aperture_radius);
    }

    activate_elens(opts, stepper_sptr, rank);

    stepper_sptr->get_lattice_simulator().update();
    if (rank == 0) {
        print_lattice(stepper_sptr->get_lattice_simulator().get_lattice_sptr(), "active_lattice.txt");
        print_beamline(stepper_sptr, "active_beamline.txt");
    }

#if 0 // tunes disabled because the lens can make the lattice unstable
    std::pair<double, double> new_tunes(stepper_sptr->get_lattice_simulator().get_both_tunes(true));
    if (rank == 0) {
        std::cout << "Tunes with electron lens active: hor: " << new_tunes.first << ", ver: " << new_tunes.second << std::endl;
    }
#endif


    Propagator propagator(stepper_sptr);
    propagator.set_checkpoint_period(opts.checkpointperiod);
    propagator.set_concurrent_io(opts.concurrentio);

    // adaptive elens setup
    Elens_actions elens_actions(opts.elensadaptive);

    propagator.propagate(bunch_train_simulator,elens_actions,opts.num_turns, opts.maxturns, opts.verbosity);

}


int
main(int argc, char **argv)
{
   MPI_Init(&argc, &argv);
   Cxx_offdiag_options opts(argc, argv);
  run(opts);
   MPI_Finalize();
    return 0;
}
