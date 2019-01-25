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

#include "ramp_actions.h"

Ramp_actions::Ramp_actions() : data(0)
{}

Ramp_actions::Ramp_actions(int initdata):
    data(initdata)
{}

Ramp_actions::~Ramp_actions()
{}

void Ramp_actions::turn_end_action(Stepper & stepper, Bunch & bunch, int turn_num)
{
    std::cout << "data: " << data << " turn_end_action, turn: " << turn_num << ", bunch num particles: " << bunch.get_total_num() << std::endl;
    ++data;
}


template<class Archive>
void Ramp_actions::serialize(Archive &ar, const unsigned int version)
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
    boost::serialization::void_cast_register<Ramp_actions, Propagate_actions>(
        static_cast<Ramp_actions *>(NULL),
        static_cast<Propagate_actions *>(NULL));

    ar & BOOST_SERIALIZATION_NVP(data);
}

template
void Ramp_actions::serialize<boost::archive::binary_oarchive >(
                                                               boost::archive::binary_oarchive &ar, const unsigned int version);

template
void Ramp_actions::serialize<boost::archive::xml_oarchive >(
                                                            boost::archive::xml_oarchive &ar, const unsigned int version);

template
void Ramp_actions::serialize<boost::archive::binary_iarchive >(
                                                               boost::archive::binary_iarchive &ar, const unsigned int version);

template
void Ramp_actions::serialize<boost::archive::xml_iarchive >(
                                                            boost::archive::xml_iarchive &ar, const unsigned int version);

BOOST_CLASS_EXPORT_IMPLEMENT(Ramp_actions)
