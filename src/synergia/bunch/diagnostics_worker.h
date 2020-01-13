#ifndef DIAGNOSTICS_WORKER_H
#define DIAGNOSTICS_WORKER_H


#include <memory>
#include <string>

#include <cereal/types/memory.hpp>

#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/diagnostics_writer.h"

class Bunch;

class Diagnostics_worker
{

private:

    std::shared_ptr<Diagnostics> diag;
    Diagnostics_writer writer;

public:

    // default constructor for serialization only
    Diagnostics_worker()
        : diag(), writer()
    { }

    // construct a diag worker with given type of diag and filename
    // specialization is provided for s_p<Diagnostics> so the python
    // interface can register
    template<class DiagCal>
    Diagnostics_worker(
            DiagCal const& diag,
            std::string const& filename,
            std::string const& local_dir = "" )
        : diag(std::make_shared<DiagCal>(diag))
        , writer(filename, local_dir)
    { }

    // for registering from python only
    Diagnostics_worker(
            std::shared_ptr<Diagnostics> const& diag,
            std::string const& filename,
            std::string const& local_dir = "" )
        : diag(diag)
        , writer(filename, local_dir)
    { }

    std::string type() const;

    void update(Bunch const& bunch);
    void write();

    void update_and_write(Bunch const& bunch)
    { update(bunch); write(); }

private:

    friend class cereal::access;

    template<class AR>
    void serialize(AR & ar)
    {
        ar(CEREAL_NVP(diag));
        ar(CEREAL_NVP(writer));
    }
};


#endif
