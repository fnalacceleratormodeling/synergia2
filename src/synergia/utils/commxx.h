#ifndef COMMXX_H_
#define COMMXX_H_

#include "mpi.h"

#include <vector>
#include <memory>

#include <cereal/types/memory.hpp>

/// Commxx is a wrapper around MPI communicator (MPI_Comm) objects.
///
/// Python:  (jfa: needs update) The equivalent functionality is provided by mpi4py Comm objects.
/// mpi4py Comm objects may be passed from python anywhere a Commxx object is
/// expected.

enum class comm_type
{ null, world, regular };

class Commxx : public std::enable_shared_from_this<Commxx>
{

public:

    static const Commxx World;
    static const Commxx Null;

    static int world_rank()
    { int r; MPI_Comm_rank(MPI_COMM_WORLD, &r); return r; }

    static int world_size()
    { int s; MPI_Comm_size(MPI_COMM_WORLD, &s); return s; }

private:

    std::shared_ptr<const MPI_Comm> comm;
    std::shared_ptr<const Commxx> parent_comm;

    comm_type type;
    int color, key;

private:

    Commxx(std::shared_ptr<const Commxx> const& parent, int color, int key);
    void construct();

public:

    /// Construct a Commxx object of MPI_COMM_WORLD or MPI_COMM_NULL
    Commxx(comm_type type = comm_type::world);

    operator MPI_Comm() const
    { if (comm) return *comm; else return MPI_COMM_NULL; }

    /// get communicator type
    comm_type get_type() const { return type; }

    /// Get communicator rank
    int get_rank() const;
    int rank() const { return get_rank(); }

    /// Get communicator size
    int get_size() const;
    int size() const { return get_size(); }

    /// Test to see if the communicator contains this rank
    bool has_this_rank() const { return (bool)comm; }
    bool is_null() const { return !(bool)comm; }

    /// is this the root Commxx object
    bool is_root() const { return type != comm_type::regular; }

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

    template<class AR>
    void save(AR & ar) const
    {
        ar(CEREAL_NVP(parent_comm));
        ar(CEREAL_NVP(type));
        ar(CEREAL_NVP(color));
        ar(CEREAL_NVP(key));
    }

    template<class AR>
    void load(AR & ar)
    { 
        ar(CEREAL_NVP(parent_comm));
        ar(CEREAL_NVP(type));
        ar(CEREAL_NVP(color));
        ar(CEREAL_NVP(key));

        switch(type)
        { 
        case comm_type::null: comm.reset(); break;
        case comm_type::world: comm.reset(new MPI_Comm(MPI_COMM_WORLD)); break;
        case comm_type::regular: construct(); break;
        }
    }
};

bool operator== (Commxx const & comm1, Commxx const & comm2);

inline bool operator!= (Commxx const & comm1, Commxx const & comm2)
{ return !(comm1 == comm2); }

#endif /* COMMXX_H_ */
