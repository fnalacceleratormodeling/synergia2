#ifndef SYNERIGA_LATTICE_LATTICE_TREE_H
#define SYNERIGA_LATTICE_LATTICE_TREE_H

#include "synergia/lattice/madx.h"

class Lattice_tree
{
public:

    explicit Lattice_tree(synergia::MadX const& madx)
        : mx(madx)
    { }

    void set_variable(std::string const& name, double val);
    void set_variable(std::string const& name, std::string const& val);

    void print() const;

//private:

    synergia::MadX mx;
};


#endif
