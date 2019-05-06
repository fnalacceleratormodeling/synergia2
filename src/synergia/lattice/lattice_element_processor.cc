
#include "lattice_element_processor.h"



Lattice_element
Lattice_element_processor::process(Lattice_element const & element)
{
    auto t = element.get_type();
    Lattice_element e(element);

    switch(element.get_type())
    {
    case element_type::drift:      process_drift(e); break;
    case element_type::quadrupole: process_quad(e); break;
    default: break;
    }

    return e;
}

void Lattice_element_processor::process_drift(Lattice_element & e)
{
    e.set_default_double_attribute("l", 0.0);
}

void Lattice_element_processor::process_quad(Lattice_element & e)
{
    e.set_default_string_attribute("propagator_type", "yoshida");
    e.set_default_double_attribute("yoshida_order", 2.0);
}
