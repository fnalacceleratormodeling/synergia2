#ifndef FF_ELEMENT_H
#define FF_ELEMENT_H

#include <beamline/JetParticle.h>
#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice_element_slice.h"

class FF_element
{
public:
    FF_element()
    : steps(default_steps), order(default_order) { };

    virtual ~FF_element();

    virtual void apply(Lattice_element_slice const& slice, JetParticle & jet_particle) = 0;
    virtual void apply(Lattice_element_slice const& slice, Bunch & bunch) = 0;

    void set_yoshida_steps(int s)
    { steps = s; }

    void set_yoshida_order(int o)
    { order = o; }

    int get_yoshida_steps() const
    { return steps; }

    int get_yoshida_order() const
    { return order; }


    // static members
    static void set_default_yoshida_steps(int s)
    { default_steps = s; }

    static void set_default_yoshida_order(int o)
    { default_order = o; }

    static int get_default_yoshdia_steps()
    { return default_steps; }

    static int get_default_yoshida_order()
    { return default_order; }

    template<class Archive>
        void serialize(Archive & ar, const unsigned int version);

protected:
    static int default_steps;
    static int default_order;

    int steps;
    int order;

    bool close_to_zero(double v)
    { return fabs(v) < 1e-13; }
};

typedef boost::shared_ptr<FF_element > FF_element_sptr; // syndoc:include

#endif // FF_ELEMENT_H
