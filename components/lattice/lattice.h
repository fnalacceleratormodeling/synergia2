#ifndef LATTICE_H_
#define LATTICE_H_

#include <string>
#include <list>


#include "components/lattice/lattice_element.h"
#include "components/foundation/reference_particle.h"

typedef std::list<Lattice_element > Lattice_element_list;

class Lattice
{
private:
    std::string name;
    Reference_particle *reference_particle_ptr;
    bool reference_particle_allocated;
    std::list<Lattice_element > elements;

public:
    Lattice(std::string const& name);
    std::string const&
    get_name() const;
    void
    set_reference_particle(Reference_particle const& reference_particle);
    bool
    has_reference_particle() const;
    Reference_particle const&
    get_reference_particle() const;
    void
    append(Lattice_element const& element);
    Lattice_element_list &
    get_elements();
    double
    get_length() const;
    double
    get_total_angle() const;
    ~Lattice();
};

#endif /* LATTICE_H_ */
