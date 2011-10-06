#include "chef_lattice_section.h"
#include <stdexcept>

Chef_lattice_section::Chef_lattice_section(Chef_lattice_sptr chef_lattice_sptr) :
    chef_lattice_sptr(chef_lattice_sptr), begin_index(no_index),
            end_index(no_index)
{
}

Chef_lattice_section::Chef_lattice_section(Chef_lattice_sptr chef_lattice_sptr,
        int begin, int end) :
    chef_lattice_sptr(chef_lattice_sptr), begin_index(begin), end_index(end)
{
}

Chef_lattice_section::Chef_lattice_section()
{
}

void
Chef_lattice_section::extend(Chef_lattice_section const& chef_lattice_section)
{
    extend(chef_lattice_section.begin_index, chef_lattice_section.end_index);
}

void
Chef_lattice_section::extend(int begin_index, int end_index)
{
    if (this->begin_index == no_index) {
        this->begin_index = begin_index;
    } else {
        if ((begin_index != this->end_index) && (begin_index != this->end_index
                + 1)) {
            std::cout << "jfa: Chef_lattice_section::extend begin_index = "
                    << begin_index << " this->end_index = " << this->end_index
                    << std::endl;
            throw std::runtime_error(
                    "Chef_lattice_section::extend: invalid begin_index");
        } else {
            this->begin_index = begin_index;
        }
    }

    if (this->end_index == no_index) {
        this->end_index = end_index;
    } else {
        if (end_index < begin_index) {
            throw std::runtime_error(
                    "Chef_lattice_section::extend: end_index less than begin_index");
        }
        this->end_index = end_index;
    }
}

int
Chef_lattice_section::get_begin_index() const
{
    return begin_index;
}

int
Chef_lattice_section::get_end_index() const
{
    return end_index;
}

Chef_lattice_section::iterator
Chef_lattice_section::begin()
{
    if (begin_index == no_index) {
        throw std::runtime_error(
                "Chef_lattice_section::begin: undefined section");
    }
    return chef_lattice_sptr->get_sliced_beamline_iterator(begin_index);
}

Chef_lattice_section::iterator
Chef_lattice_section::end()
{
    if (end_index == no_index) {
        throw std::runtime_error("Chef_lattice_section::end: undefined section");
    }
    return chef_lattice_sptr->get_sliced_beamline_iterator(end_index);
}

Chef_lattice_section::const_iterator
Chef_lattice_section::begin() const
{
    if (begin_index == no_index) {
        throw std::runtime_error(
                "Chef_lattice_section::begin: undefined section");
    }
    return chef_lattice_sptr->get_sliced_beamline_const_iterator(begin_index);
}

Chef_lattice_section::const_iterator
Chef_lattice_section::end() const
{
    if (end_index == no_index) {
        throw std::runtime_error("Chef_lattice_section::end: undefined section");
    }
    return chef_lattice_sptr->get_sliced_beamline_const_iterator(end_index);
}

bool
Chef_lattice_section::empty() const
{
    return (begin_index == no_index);
}

void
Chef_lattice_section::clear()
{
    begin_index = no_index;
    end_index = no_index;
}
