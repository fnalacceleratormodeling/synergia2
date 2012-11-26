#ifndef GENERALIZED_DIAGNOSTICS_H_
#define GENERALIZED_DIAGNOSTICS_H_

#include <string>
#include "boost/shared_ptr.hpp"
#include "synergia/utils/serialization.h"

/// Generalized_diagnostics is an abstract base class for bunch and lattice
/// diagnostics classes
class Generalized_diagnostics
{
private:
    std::string name;
public:
    Generalized_diagnostics(std::string const& name);
    // Default constructor for serialization use only
    Generalized_diagnostics();
    /// return diagnostics type
    std::string const &
    get_name() const;
    /// Multiple serial diagnostics can be written to a single file.
    virtual bool
    is_serial() const = 0;
    /// Update the diagnostics
    virtual void
    update() = 0;
    /// Write the diagnostics to the file
    virtual void
    write() = 0;
    /// Update the diagnostics and write them to the file
    virtual void
    update_and_write();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(name);
        }
    virtual
    ~Generalized_diagnostics()
    {
    }
    ;
};

typedef boost::shared_ptr<Generalized_diagnostics >
        Generalized_diagnostics_sptr; // syndoc:include

#endif /* GENERALIZED_DIAGNOSTICS_H_ */
