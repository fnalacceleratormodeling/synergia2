#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <string>
#include <list>
#include <boost/shared_ptr.hpp>

#include "synergia/bunch/bunch.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/simulation/independent_operation.h"
#include "synergia/simulation/operation_extractor.h"
#include "synergia/simulation/aperture_operation_extractor.h"
#include "synergia/foundation/multi_diagnostics.h"
#include "synergia/bunch/train.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/logger.h"

class Step;

class Operator
{
private:
    std::string name, type;
public:
    Operator(std::string const& name, std::string const& type);
    /// Default constructor for serialization use only
    Operator();
    std::string const&
    get_name() const;
    std::string const&
    get_type() const;
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger) = 0;
#if 0
    virtual void
    apply_train(Bunch_with_diagnostics_train & bunch_diag_train,
            double time_step, Step & step);
#endif
    virtual void
    print() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Operator();
};
typedef boost::shared_ptr<Operator > Operator_sptr; // syndoc:include
typedef std::list<Operator_sptr > Operators; // syndoc:include

class Collective_operator : public Operator
{
public:
    Collective_operator(std::string const& name);
    /// Default constructor for serialization use only
    Collective_operator();
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger) = 0;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    static const char type_name[];
    virtual
    ~Collective_operator();
};
BOOST_CLASS_EXPORT_KEY(Collective_operator);
typedef boost::shared_ptr<Collective_operator > Collective_operator_sptr; // syndoc:include
typedef std::list<Collective_operator_sptr > Collective_operators; // syndoc:include

class Dummy_collective_operator : public Collective_operator
{
public:
    Dummy_collective_operator(std::string const& name);
    /// Default constructor for serialization use only
    Dummy_collective_operator();
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity, Logger & logger);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Dummy_collective_operator();
};
BOOST_CLASS_EXPORT_KEY(Dummy_collective_operator);
typedef boost::shared_ptr<Dummy_collective_operator >
        Dummy_collective_operator_sptr; // syndoc:include
typedef std::list<Dummy_collective_operator_sptr > Dummy_collective_operators; // syndoc:include

class Independent_operator : public Operator
{
private:
    Lattice_element_slices slices;
    Independent_operations operations;
    std::list<long int > operations_revisions;
    Reference_particle operations_reference_particle;
    Operation_extractor_map_sptr operation_extractor_map_sptr;
    Aperture_operation_extractor_map_sptr aperture_operation_extractor_map_sptr;
    bool have_operations;
public:
    Independent_operator(
            std::string const& name,
            Operation_extractor_map_sptr operation_extractor_map_sptr,
            Aperture_operation_extractor_map_sptr aperture_operation_extractor_map_sptr);
    /// Default constructor for serialization use only
    Independent_operator();
    void
    append_slice(Lattice_element_slice_sptr slice_sptr);
    Lattice_element_slices const&
    get_slices() const;
    void
    update_operations(Reference_particle const& reference_particle);
    bool
    need_update(Reference_particle const& reference_particle, int verbosity,
            Logger & logger);
    Independent_operations const&
    get_operations() const;
    Independent_operations &
    get_operations();
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger);
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger, Multi_diagnostics & diagnostics);
    virtual void
    print() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    static const char type_name[];
    virtual
    ~Independent_operator();
};
BOOST_CLASS_EXPORT_KEY(Independent_operator);
typedef boost::shared_ptr<Independent_operator > Independent_operator_sptr; // syndoc:include
typedef std::list<Independent_operator_sptr > Independent_operators; // syndoc:include

#endif /* OPERATOR_H_ */
