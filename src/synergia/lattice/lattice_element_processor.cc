
#include "lattice_element_processor.h"



Lattice_element
Lattice_element_processor::process(Lattice_element const & element)
{
    auto t = element.get_type();
    Lattice_element e(element);

    switch(element.get_type())
    {
    case element_type::drift:      drift(e); break;
    case element_type::sbend:      sbend(e); break;
    case element_type::quadrupole: quadrupole(e); break;
    case element_type::multipole:  multipole(e); break;
    case element_type::rfcavity:   rfcavity(e); break;
    default: break;
    }

    return e;
}

void Lattice_element_processor::drift(Lattice_element & e)
{
    e.set_default_double_attribute("l", 0.0);
}

void Lattice_element_processor::sbend(Lattice_element & e)
{
    e.set_default_string_attribute("propagator_type", "yoshida");
    e.set_default_double_attribute("l", 0.0);
    e.set_default_double_attribute("angle", 0.0);
    e.set_default_double_attribute("tilt", 0.0);
    e.set_default_double_attribute("k1", 0.0);
    e.set_default_double_attribute("e1", 0.0);
    e.set_default_double_attribute("e2", 0.0);
    e.set_default_double_attribute("fint", 0.0);
    e.set_default_double_attribute("fintx", 0.0);
    e.set_default_double_attribute("hgap", 0.0);
    e.set_default_double_attribute("k2", 0.0);
    e.set_default_double_attribute("h1", 0.0);
    e.set_default_double_attribute("h2", 0.0);
    e.set_default_double_attribute("kicks", 40.0);

    // possible higher order multipole components
    e.set_default_double_attribute("kl", 0.0); // base strength/B-rho
    e.set_default_double_attribute("a1", 0.0); // skew quad
    e.set_default_double_attribute("a2", 0.0); // skew sextupole
    e.set_default_double_attribute("a3", 0.0); // skew octupole
    e.set_default_double_attribute("a4", 0.0); // skew decapole
    e.set_default_double_attribute("a5", 0.0); // skew dodecapole
    e.set_default_double_attribute("a6", 0.0); // skew tetradecapole
    e.set_default_double_attribute("a7", 0.0); // skew hexdecapole
    e.set_default_double_attribute("b1", 0.0); // quad
    e.set_default_double_attribute("b2", 0.0); // sextupole
    e.set_default_double_attribute("b3", 0.0); // octopole
    e.set_default_double_attribute("b4", 0.0); // decapole
    e.set_default_double_attribute("b5", 0.0); // dodecapole
    e.set_default_double_attribute("b6", 0.0); // tetradecapole
    e.set_default_double_attribute("b7", 0.0); // hexdecapole
    e.set_default_double_attribute("entry_edge_kick", 1.0);
    e.set_default_double_attribute("exit_edge_kick", 1.0);
}

void Lattice_element_processor::quadrupole(Lattice_element & e)
{
    e.set_default_string_attribute("propagator_type", "yoshida");
    e.set_default_double_attribute("yoshida_order", 2.0);
}

void Lattice_element_processor::multipole(Lattice_element & e)
{
    e.set_default_double_attribute("tilt", 0.0);
    //e.set_default_vector_attribute("knl", std::vector<double >(0));
    //e.set_default_vector_attribute("ksl", std::vector<double >(0));
}

void Lattice_element_processor::rfcavity(Lattice_element & e)
{
    e.set_default_double_attribute("l", 0.0);
    e.set_default_double_attribute("volt", 0.0);
    e.set_default_double_attribute("lag", 0.0);
    e.set_default_double_attribute("harmon", 0.0);
    e.set_default_double_attribute("shunt", 0.0);
}




