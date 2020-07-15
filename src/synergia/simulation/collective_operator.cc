#include "synergia/simulation/collective_operator.h"


const char Collective_operator::type_name[] = "collective";

Collective_operator::Collective_operator(std::string const& name) :
        Operator(name, type_name)
{
}

Collective_operator::Collective_operator()
{
}

void
Collective_operator::apply(Bunch & bunch, double time_step, Step & step, int verbosity,
        Diagnosticss const& per_operation_diagnosticss, Logger & logger)
{
    apply(bunch, time_step, step, verbosity, logger);
}

template<class Archive>
    void
    Collective_operator::serialize(Archive & ar, const unsigned int version)
    {
        ar &
        BOOST_SERIALIZATION_BASE_OBJECT_NVP(Operator);
    }

template
void
Collective_operator::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Collective_operator::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Collective_operator::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Collective_operator::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Collective_operator::~Collective_operator()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT (Collective_operator)

Dummy_collective_operator::Dummy_collective_operator(std::string const& name) :
        Collective_operator(name)
{
}

Dummy_collective_operator::Dummy_collective_operator()
{
}

Dummy_collective_operator *
Dummy_collective_operator::clone()
{
    return new Dummy_collective_operator(*this);
}

void
Dummy_collective_operator::apply(Bunch & bunch, double time_step, Step & step,
        int verbosity, Logger & logger)
{
}

template<class Archive>
    void
    Dummy_collective_operator::serialize(Archive & ar,
            const unsigned int version)
    {
        ar &
        BOOST_SERIALIZATION_BASE_OBJECT_NVP(Collective_operator);
    }

template
void
Dummy_collective_operator::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Dummy_collective_operator::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Dummy_collective_operator::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Dummy_collective_operator::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Dummy_collective_operator::~Dummy_collective_operator()
{
}
BOOST_CLASS_EXPORT_IMPLEMENT (Dummy_collective_operator)

