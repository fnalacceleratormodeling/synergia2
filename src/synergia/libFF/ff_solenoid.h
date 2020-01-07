#ifndef FF_SOLENOID_H
#define FF_SOLENOID_H

#include "synergia/libFF/ff_algorithm.h"
#include "synergia/libFF/ff_element.h"

class FF_solenoid : public FF_element
{
private:
    double get_reference_cdt_drift(double length, Reference_particle & reference_particle);
    double get_reference_cdt_solenoid(double length, Reference_particle & reference_particle,
            bool in_edge, bool out_edge, double ks, double kse, double ksl);
public:
    FF_solenoid();
    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle);
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch);
    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);
    virtual ~FF_solenoid();
};

typedef boost::shared_ptr<FF_solenoid > FF_solenoid_sptr;

#endif // FF_SOLENOID_H