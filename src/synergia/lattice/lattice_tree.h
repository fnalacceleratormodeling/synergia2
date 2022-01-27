#ifndef SYNERIGA_LATTICE_LATTICE_TREE_H
#define SYNERIGA_LATTICE_LATTICE_TREE_H

#include "synergia/lattice/madx.h"
#include "synergia/lattice/mx_parse.h"

#include "synergia/utils/cereal.h"

class Lattice_tree
{
public:

    // default ctor for serialization
    Lattice_tree() : mx()
    { }

    explicit Lattice_tree(synergia::MadX const& madx)
        : mx(madx)
    { }

    void set_variable(std::string const& name, double val);
    void set_variable(std::string const& name, std::string const& val);

    void print() const;

public:

    synergia::MadX mx;

private:

    friend class cereal::access;


    template<class Archive>
    void save(Archive & ar) const
    {
        std::string madx = mx.to_madx();
        ar(madx);
    }

    template<class Archive>
    void load(Archive & ar)
    {
        std::string madx;
        ar(madx);

        parse_madx(madx, mx);
    }
};


#endif
