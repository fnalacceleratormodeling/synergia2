#ifndef COLLECTIVE_OPERATOR_H_
#define COLLECTIVE_OPERATOR_H_

#include "synergia/simulation/operator.h"

class Collective_operator : public Operator
{
public:
    Collective_operator(std::string const& name);
    /// Default constructor for serialization use only
    Collective_operator();
    virtual Collective_operator *
    clone() = 0;
    using Operator::apply;
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Diagnosticss const& per_operation_diagnosticss, Logger & logger);
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
    virtual Dummy_collective_operator *
    clone();
    using Collective_operator::apply;
    virtual void
    apply(Bunch & bunch, double time_step, Step & step, int verbosity,
            Logger & logger);
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



#endif
