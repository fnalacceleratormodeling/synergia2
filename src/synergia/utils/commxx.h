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

    static const Commxx world;

private:

    std::shared_ptr<MPI_Comm> comm;
    std::weak_ptr<MPI_Comm> parent_comm;

private:

    Commxx(std::shared_ptr<MPI_Comm> c);

public:

    /// Construct a Commxx object using MPI_COMM_WORLD
    Commxx();
    Commxx(MPI_Comm const & comm, comm_create_kind kind);

    operator MPI_Comm() const
    { if (comm) return *comm; else return MPI_COMM_NULL; }

    /// Get communicator rank
    int get_rank() const;
    int rank() const { return get_rank(); }

    /// Get communicator size
    int get_size() const;
    int size() const { return get_size(); }

    /// Test to see if the communicator contains this rank
    bool has_this_rank() const { return (bool)comm; }
    bool is_null() const { return !(bool)comm; }

    // get the parent communicator if available
    Commxx parent() const;

    // the rank in the context of parent communicator
    // int rank_in_parent_comm() const;

    // the ranks (in the context of parent communicator) of 
    // all processors in the current communicator
    // std::vector<int> ranks_in_parent_comm() const;

    Commxx dup() const;
    Commxx split(int color) const;
    Commxx split(int color, int key) const;
    Commxx divide(int subgroup_size) const;
    Commxx group(std::vector<int> const & ranks) const;
};

bool operator== (Commxx const & comm1, Commxx const & comm2);

inline bool operator!= (Commxx const & comm1, Commxx const & comm2)
{ return !(comm1 == comm2); }

#endif /* COMMXX_H_ */
