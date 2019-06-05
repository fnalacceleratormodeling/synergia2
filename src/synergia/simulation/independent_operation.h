#ifndef INDEPENDENT_OPERATION_H_
#define INDEPENDENT_OPERATION_H_

//#include "synergia/simulation/fast_mapping.h"
//#include "synergia/simulation/chef_propagator.h"

#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/libFF/ff_element.h"

#include "synergia/utils/logger.h"



class Independent_operation
{
private:

    std::string type;

private:

    virtual void print_impl(Logger & logger) const = 0;

public:

    Independent_operation(std::string const & type)
    : type(type)
    { }

    std::string const & get_type() const
    { return type; }

    virtual void apply(Bunch & bunch, Logger & logger) 
    { }

    void print(Logger & logger) const
    { 
        logger(LoggerV::DEBUG) << "\ttype = " << type << ", ";
        print_impl(logger); 
    }
};


class LibFF_operation : public Independent_operation
{
private:

    std::vector<std::pair<std::unique_ptr<FF_element>, Lattice_element_slice>>
        libff_element_slices;

private:

    void print_impl(Logger & logger) const override
    { }

public:

    LibFF_operation(std::vector<Lattice_element_slice> const & slices);

    virtual void apply(Bunch & bunch, Logger & logger);
};


#if 0
class Independent_operation
{
private:
    std::string type;
public:
    Independent_operation(std::string const& type);
    /// Default constructor for serialization use only
    Independent_operation();
    std::string const&
    get_type() const;
    virtual void
    apply(Bunch & bunch, int verbosity, Logger & logger) = 0;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Independent_operation();
};

typedef boost::shared_ptr<Independent_operation > Independent_operation_sptr; // syndoc:include
typedef std::list<Independent_operation_sptr > Independent_operations; // syndoc:include

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
    apply(Bunch & bunch, int verbosity, Logger & logger);
    Fast_mapping const&
    get_fast_mapping() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Fast_mapping_operation();
};
BOOST_CLASS_EXPORT_KEY(Fast_mapping_operation);
typedef boost::shared_ptr<Fast_mapping_operation > Fast_mapping_operation_sptr; // syndoc:include

const char chef_propagate_type_name[] = "chef_propagate";
class Chef_propagate_operation : public Independent_operation
{
private:
    Chef_propagator chef_propagator;

public:
    Chef_propagate_operation(Chef_lattice_section const& chef_lattice_section);
    /// Default constructor for serialization use only
    Chef_propagate_operation();
    virtual void
    apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Chef_propagate_operation();
};
BOOST_CLASS_EXPORT_KEY(Chef_propagate_operation);
typedef boost::shared_ptr<Chef_propagate_operation >
        Chef_propagate_operation_sptr; // syndoc:include

const char libFF_type_name[] = "LibFF";
class LibFF_operation : public Independent_operation
{
private:
    Lattice_element_slices lattice_element_slices;
public:
    LibFF_operation(Lattice_element_slices const& lattice_element_slices);
    /// Default constructor for serialization use only
    LibFF_operation();
    virtual void
    apply(Bunch & bunch, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~LibFF_operation();
};
BOOST_CLASS_EXPORT_KEY(LibFF_operation);
typedef boost::shared_ptr<LibFF_operation >
        LibFF_operation_sptr; // syndoc:include
#endif

#endif /* INDEPENDENT_OPERATION_H_ */
