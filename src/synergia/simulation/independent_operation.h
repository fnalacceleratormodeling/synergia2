#ifndef INDEPENDENT_OPERATION_H_
#define INDEPENDENT_OPERATION_H_

#include <list>
#include <string>
#include <map>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/shared_ptr.hpp>

#include "synergia/simulation/fast_mapping.h"
#include "synergia/simulation/chef_propagator.h"
#include "synergia/lattice/chef_lattice.h"

class Independent_operation
{
private:
    std::string type;
public:
    Independent_operation(std::string const& type);
    /// Default constructor for serialization use only
    Independent_operation();
    std::string
    get_type() const;
    virtual void
    apply(Bunch & bunch) = 0;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(type);
        }
    virtual
    ~Independent_operation();
};

typedef boost::shared_ptr<Independent_operation > Independent_operation_sptr;
typedef std::list<Independent_operation_sptr > Independent_operations;

const char fast_mapping_type_name[] = "fast_mapping";
class Fast_mapping_operation : public Independent_operation
{
private:
    Fast_mapping mapping;

public:
    Fast_mapping_operation(Fast_mapping const& mapping);
    /// Default constructor for serialization use only
    Fast_mapping_operation();
    virtual void
    apply(Bunch & bunch);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Independent_operation);
            ar & BOOST_SERIALIZATION_NVP(mapping);
        }
    virtual
    ~Fast_mapping_operation();
};

typedef boost::shared_ptr<Fast_mapping_operation > Fast_mapping_operation_sptr;

const char chef_propagate_type_name[] = "chef_propagate";
class Chef_propagate_operation : public Independent_operation
{
private:
    Chef_propagator chef_propagator;

public:
    Chef_propagate_operation(Chef_lattice_section_sptr chef_lattice_section_sptr);
    /// Default constructor for serialization use only
    Chef_propagate_operation();
    virtual void
    apply(Bunch & bunch);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP(Independent_operation);
            ar & BOOST_SERIALIZATION_NVP(chef_propagator);
        }
    virtual
    ~Chef_propagate_operation();
};

typedef boost::shared_ptr<Chef_propagate_operation >
        Chef_propagate_operation_sptr;

#endif /* INDEPENDENT_OPERATION_H_ */
