#include "synergia/utils/catch.hpp"


#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element.h"

#include "synergia/lattice/madx_reader.h"
#include "synergia/utils/commxx.h"
#include "synergia/utils/logger.h"
#include "synergia/utils/simple_timer.h"
#include "synergia/utils/utils.h"
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/independent_stepper_elements.h"

#include <iomanip>
#include <iostream>
#include <string>


Lattice
get_lattice()
{
    static std::string fodo_madx(R"foo(
beam, particle=proton,pc=3.0;

o: drift, l=8.0;
f: quadrupole, l=2.0, k1=0.071428571428571425, a1=1.0;
d: quadrupole, l=2.0, k1=-0.071428571428571425;

fodo: sequence, l=20.0, refer=entry;
fodo_1: f, at=0.0;
fodo_2: o, at=2.0;
fodo_3: d, at=10.0;
fodo_4: o, at=12.0;
endsequence;
)foo");

    MadX_reader reader;
    reader.parse(fodo_madx);
    return reader.get_lattice("fodo");
}

Propagator
create_propagator(Lattice lattice)
{
    Independent_stepper_elements stepper(1);
    Propagator prop(lattice, stepper);
    return prop;
}

TEST_CASE("create_propagator")
{
    Lattice lattice(get_lattice());
    Propagator p(create_propagator(lattice));
}

TEST_CASE("test_prop_get_lattice_elements")
{
    Lattice lattice(get_lattice());
    Propagator p(create_propagator(lattice));
   
    int nlattice_elem = lattice.get_elements().size();
    int nprop_elem = p.get_lattice().get_elements().size();
    CHECK (nlattice_elem == nprop_elem);
}

TEST_CASE("test_modify_lattice_elements")
{
    Lattice lattice(get_lattice());
    Propagator p(create_propagator(lattice));

    std::list<Lattice_element>& elems(p.get_lattice().get_elements());
    std::list<Lattice_element>::iterator first = elems.begin();
    int orig_a1 = (*first).get_double_attribute("a1");

    int new_a1 = orig_a1 + 100.0;
    (*first).set_double_attribute("a1", new_a1);

    std::list<Lattice_element> const& newelems(p.get_lattice().get_elements());
    std::list<Lattice_element>::const_iterator newfirst = newelems.begin();
    CHECK ( (*newfirst).get_double_attribute("a1") == Approx(new_a1) );

    Propagator::Lattice_element_slices les = p.get_lattice_element_slices();
    // with Independent_stepper_elements, first slice is first element
    // check a1 attribute
    Propagator::Slice_iterator sit = les.begin();

    double test_a1 = (*sit).get_lattice_element().get_double_attribute("a1");
    CHECK (test_a1 == Approx(new_a1));
}

TEST_CASE("test_modify_reference_particle_energy")
{
    Lattice lattice(get_lattice());
    Propagator p(create_propagator(lattice));

    double orig_energy = lattice.get_reference_particle().get_total_energy();
    int new_energy = orig_energy + 1.0;
    p.get_lattice().get_reference_particle().set_total_energy(new_energy);

    CHECK ( p.get_lattice().get_reference_particle().get_total_energy() == Approx(new_energy) );

}
