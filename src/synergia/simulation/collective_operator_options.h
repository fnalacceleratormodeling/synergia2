#ifndef COLLECTIVE_OPERATOR_OPTIONS_H
#define COLLECTIVE_OPERATOR_OPTIONS_H

#include <cereal/types/polymorphic.hpp>

class Collective_operator;

struct CO_options
{
    virtual CO_options * clone() const = 0;
    virtual Collective_operator * create_operator() const = 0;

    template<class Archive>
    void serialize(Archive & ar) { }
};

#include <cereal/archives/json.hpp>



#endif
