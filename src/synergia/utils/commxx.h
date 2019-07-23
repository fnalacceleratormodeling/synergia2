#ifndef COMMXX_H_
#define COMMXX_H_

#include <vector>
#include <stdexcept>
#include <memory>
#include "mpi.h"

//#include "synergia/utils/cereal.h"

enum class comm_create_kind { duplicate, take_ownership, attach };

/// Commxx is a wrapper around MPI communicator (MPI_Comm) objects.
///
/// Python:  (jfa: needs update) The equivalent functionality is provided by mpi4py Comm objects.
/// mpi4py Comm objects may be passed from python anywhere a Commxx object is
/// expected.

class Commxx
{
public:

    static Commxx Commxx_world;

private:

    std::shared_ptr<MPI_Comm> comm;
    std::weak_ptr<MPI_Comm> parent_comm;
    //Commxx const & parent;

    //bool has_this_rank_;

public:

    /// Construct a Commxx object using MPI_COMM_WORLD
    Commxx();
    Commxx(MPI_Comm const & comm, comm_create_kind kind);

    operator MPI_Comm() const
    { if (comm) return *comm; else return MPI_COMM_NULL; }

    //Commxx const & get_parent() const
    //{ return parent; }

    /// Get communicator rank
    int get_rank() const;

    /// Get communicator size
    int get_size() const;

    /// Test to see if the communicator contains this rank
    bool has_this_rank() const { return (bool)comm; }

    /// Extract the MPI_comm object wrapped by the Commxx instance.
    //MPI_Comm get() const;
};

#endif /* COMMXX_H_ */
