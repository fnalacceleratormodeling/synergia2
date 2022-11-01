#include <stdexcept>

#include "commxx.h"

namespace {
    struct comm_free {
        void
        operator()(MPI_Comm* comm)
        {
            if (comm == nullptr) return;
            if (*comm == MPI_COMM_NULL) return;

            MPI_Comm_free(comm);
            delete comm;
        }
    };

}

const Commxx Commxx::World(comm_type::world);
const Commxx Commxx::Null(comm_type::null);

Commxx::Commxx(comm_type type)
    : comm(type == comm_type::null ? nullptr : new MPI_Comm(MPI_COMM_WORLD))
    , parent_comm()
    , type(type)
    , color(0)
    , key(0)
{
    if (type == comm_type::regular)
        throw std::runtime_error("only null or mpi_comm_world can be created");
}

void
Commxx::construct()
{
    if (!parent_comm || parent_comm->is_null())
        throw std::runtime_error(
            "invalid parent communicator while in construct");

    MPI_Comm newcomm;
    MPI_Comm_split(*(parent_comm->comm), color, key, &newcomm);

    // do not construct the null communicator
    // always take the ownership
    if (newcomm != MPI_COMM_NULL)
        comm.reset(new MPI_Comm(newcomm), comm_free());
}

Commxx::Commxx(std::shared_ptr<const Commxx> const& parent, int color, int key)
    : comm()
    , parent_comm(std::move(parent))
    , type(comm_type::regular)
    , color(color)
    , key(key)
{
    construct();
}

int
Commxx::get_rank() const
{
    if (type == comm_type::null)
        throw std::runtime_error("Cannot get_rank() for a null commxx");

    int error, rank;
    error = MPI_Comm_rank(*comm, &rank);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_rank");
    }
    return rank;
}

int
Commxx::get_size() const
{
    if (type == comm_type::null)
        throw std::runtime_error("Cannot get_size() for a null commxx");

    int error, size;
    error = MPI_Comm_size(*comm, &size);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_size");
    }
    return size;
}

Commxx
Commxx::parent() const
{
    if (!parent_comm) throw std::runtime_error("Invalid parent communicator");

    return *parent_comm;
}

Commxx
Commxx::dup() const
{
    if (is_null()) throw std::runtime_error("dup from a null comm");
    return Commxx(shared_from_this(), 0, rank());
}

Commxx
Commxx::split(int color) const
{
    if (is_null()) throw std::runtime_error("split from a null comm");
    return Commxx(shared_from_this(), color, rank());
}

Commxx
Commxx::split(int color, int key) const
{
    if (is_null()) throw std::runtime_error("split from a null comm");
    return Commxx(shared_from_this(), color, key);
}

Commxx
Commxx::divide(int subgroup_size) const
{
    if (is_null()) throw std::runtime_error("divide from a null comm");

    if (size() < subgroup_size) return dup();
    if (size() % subgroup_size != 0) {
        throw std::runtime_error(
            "Commxx::divide(): size must be divisible by subgroup_size");
    }

    int color = rank() / subgroup_size;
    return split(color);
}

Commxx
Commxx::group(std::vector<int> const& ranks) const
{
    if (is_null()) throw std::runtime_error("group from a null comm");

    int r = rank();
    int color = (std::find(ranks.begin(), ranks.end(), r) != ranks.end()) ?
                    0 :
                    MPI_UNDEFINED;

    return Commxx(shared_from_this(), color, r);
}

bool
operator==(Commxx const& comm1, Commxx const& comm2)
{
    // both null
    if (comm1.is_null() && comm2.is_null()) return true;

    // one is null, the other is not
    if (comm1.is_null() || comm2.is_null()) return false;

    // both not null
    int result;
    MPI_Comm_compare((MPI_Comm)comm1, (MPI_Comm)comm2, &result);
    return result == MPI_IDENT;
}
