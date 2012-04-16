#ifndef CHEF_LATTICE_SECTION_H_
#define CHEF_LATTICE_SECTION_H_

#include "synergia/lattice/chef_lattice.h"
#include "synergia/lattice/chef_lattice_section_fwd.h"

class Chef_lattice_section
{
private:
    int begin_index, end_index;
    Chef_lattice_sptr chef_lattice_sptr;
    static const int no_index = -999;

public:
    typedef beamline::iterator iterator;
    typedef beamline::const_iterator const_iterator;

    Chef_lattice_section(Chef_lattice_sptr chef_lattice_sptr);
    Chef_lattice_section(Chef_lattice_sptr chef_lattice_sptr, int begin_index,
            int end_index);
    /// Default constructor for serialization use only
    Chef_lattice_section();
    void
    extend(Chef_lattice_section const& chef_lattice_section);
    void
    extend(int begin_index, int end_index);
    int
    get_begin_index() const;
    int
    get_end_index() const;
    iterator
    begin();
    iterator
    end();
    const_iterator
    begin() const;
    const_iterator
    end() const;
    bool
    empty() const;
    void
    clear();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};

#endif /* CHEF_LATTICE_SECTION_H_ */
