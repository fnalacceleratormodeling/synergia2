#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/serialization.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/simulation/propagate_actions.h"

#define DEBUG 0

#include "elens_actions.h"

Elens_actions::Elens_actions() : adaptive(false)
{}

Elens_actions::Elens_actions(bool elensadaptive):
    adaptive(elensadaptive)
{
#if DEBUG
    std::cout << "elens_actions: elens_actions created, adaptive is " <<
        (adaptive?"ON":"OFF") << std::endl;
#endif
}

Elens_actions::~Elens_actions()
{
}

void Elens_actions::turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num)
{

#if DEBUG
    if (bunch.get_comm().get_rank() == 0) {
        std::cout << "elens_actions:turn_end_action turn: " << turn_num << ", adaptive: " << adaptive << std::endl;
    }
#endif

    if (adaptive) {
        set_elens_adaptive_radius(stepper, bunch);
    }
}

void Elens_actions::turn_end_action(Stepper & stepper, Bunch_train & bunch_train, int turn_num)
{

#if DEBUG
    if (bunch_train.get_parent_comm_sptr()->get_rank() == 0) {
        std::cout << "elens_actions:turn_end_action turn: " << turn_num << ", adaptive: " << adaptive << std::endl;
    }
#endif

    if (adaptive) {
        set_elens_adaptive_radius(stepper, bunch_train);
    }
}

void Elens_actions::first_action(Stepper & stepper, Bunch & bunch)
{

#if DEBUG
    if (bunch.get_comm().get_rank() == 0) {
        std::cout << "elens_actions:first_action: adaptive: " << adaptive << std::endl;
    }
#endif

    if (adaptive) {
        set_elens_adaptive_radius(stepper, bunch);
    }
}

void Elens_actions::first_action(Stepper & stepper, Bunch_train & bunch_train)
{

#if DEBUG
    if (bunch_train.get_parent_comm_sptr()->get_rank() == 0) {
        std::cout << "elens_actions:first_action turn: adaptive: " << adaptive << std::endl;
    }
#endif

    if (adaptive) {
        set_elens_adaptive_radius(stepper, bunch_train);
    }
}


void
Elens_actions::set_elens_adaptive_radius(Stepper &stepper, Bunch_train &bunch_train)
{
    // only using the first bunch of the train
    Bunch & bunch(*bunch_train.get_bunches()[0]);
    set_elens_adaptive_radius(stepper, bunch);
}


void
Elens_actions::set_elens_adaptive_radius(Stepper &stepper, Bunch &bunch)
{
    // if adaptive:
    //     calculate the bunch RMS
    //     set all the elens radii to that value
    //     update the lattice simulator

  
    MArray1d bunch_mean(Core_diagnostics::calculate_spatial_mean(bunch));
    MArray1d bunch_std(Core_diagnostics::calculate_spatial_std(bunch, bunch_mean));

    // use geometric mean of x and y
    double elensradius = std::sqrt(bunch_std[0]*bunch_std[1]);
#if DEBUG
    if (bunch.get_comm().get_rank() == 0) {
        std::cout << "elens_actions: set_elens_adaptive_radius: elens radius: " << elensradius << std::endl;
    }
#endif

    int elenscnt = 0;
    Lattice_sptr lattice_sptr(stepper.get_lattice_simulator().get_lattice_sptr());
    for (Lattice_elements::iterator leit=lattice_sptr->get_elements().begin();
         leit!=lattice_sptr->get_elements().end(); ++leit) {
        if ((*leit)->get_type() == "elens") {
            (*leit)->set_double_attribute("radius", elensradius);
            ++elenscnt;
        }
    }
#if DEBUG
    if (bunch.get_comm().get_rank() == 0) {
        std::cout << "elens_actions: set_elens_adaptive_radius: " << elenscnt << " elens set"<< std::endl;
    }
#endif

    if (elenscnt > 0) {
        stepper.get_lattice_simulator().update();
#if DEBUG
        if (bunch.get_comm().get_rank() == 0) {
            std::cout << "elens_actions: set_elens_adaptive_radius: updating lattice_simulator" << std::endl;
        }
#endif
    }
}


template<class Archive>
void Elens_actions::serialize(Archive &ar, const unsigned int version)
{

    // The next line builds fine but doesn't appear to be necessary.  That's
    // probably become the base class has no members so there is nothing
    // to save.
    ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Propagate_actions);

    // In order to register the relationship between the derived and class
    // classes, I have to register it.  Since the base class doesn't have
    // its own serialization method, I have to do it explicitly.  Recipe is
    // from file boost/libs/serialization/doc/serialization.html#runtimecasting
    // method 2 : explicitly register base/derived relationship
    boost::serialization::void_cast_register<Elens_actions, Propagate_actions>(
        static_cast<Elens_actions *>(NULL),
        static_cast<Propagate_actions *>(NULL));

    ar & BOOST_SERIALIZATION_NVP(adaptive);
}

template
void Elens_actions::serialize<boost::archive::binary_oarchive >(
                                                               boost::archive::binary_oarchive &ar, const unsigned int version);

template
void Elens_actions::serialize<boost::archive::xml_oarchive >(
                                                            boost::archive::xml_oarchive &ar, const unsigned int version);

template
void Elens_actions::serialize<boost::archive::binary_iarchive >(
                                                               boost::archive::binary_iarchive &ar, const unsigned int version);

template
void Elens_actions::serialize<boost::archive::xml_iarchive >(
                                                            boost::archive::xml_iarchive &ar, const unsigned int version);

BOOST_CLASS_EXPORT_IMPLEMENT(Elens_actions)
