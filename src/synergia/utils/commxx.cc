#include "commxx.h"
#include <stdexcept>
#include <climits>
#include <algorithm>


namespace
{
    struct comm_free
    {
        void operator() (MPI_Comm * comm)
        {
            if ( comm == nullptr) return;
            if (*comm == MPI_COMM_NULL) return;

            MPI_Comm_free(comm);
            delete comm;
        }
    };

    size_t hash(const char * s)
    {
        size_t h = 37062913;
        while (*s) h = h * 101 + (unsigned char) *s++;
        return h;
    }
}



#if 0
void
Commxx::construct(MPI_Comm const& parent_mpi_comm)
{
    MPI_Comm temp_comm;
    int error;

    if (ranks.size() > 0) {
        MPI_Group parent_group, group;
        error = MPI_Comm_group(parent_mpi_comm, &parent_group);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in Commxx(MPI_Comm_group)");
        }
        error = MPI_Group_incl(parent_group, ranks.size(), &ranks[0], &group);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in Commxx(MPI_Group_incl)");
        }
        error = MPI_Comm_create(parent_mpi_comm, group, &temp_comm);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in Commxx(MPI_Comm_create)");
        }
        error = MPI_Group_free(&parent_group);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error(
                    "MPI error in Commxx(MPI_Group_free(parent))");
        }
        error = MPI_Group_free(&group);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in Commxx(MPI_Group_free)");
        }

        has_this_rank_ = false;
        const int no_rank = -1;
        int this_rank = no_rank;
        if (!parent_sptr) {
            this_rank = Commxx().get_rank();
        } else {
            if (parent_sptr->has_this_rank()) {
                this_rank = parent_sptr->get_rank();
            }
        }
        if (this_rank != no_rank) {
            for (std::vector<int >::const_iterator it = ranks.begin(); it
                    != ranks.end(); ++it) {
                if ((*it) == this_rank) {
                    has_this_rank_ = true;
                }
            }
        }

    } else {
        temp_comm = parent_mpi_comm;
        has_this_rank_ = true;
    }

    if (per_host && has_this_rank_) {
        char name[MPI_MAX_PROCESSOR_NAME];
        int name_len;
        MPI_Get_processor_name(name, &name_len);

        int color = hash(name) % INT_MAX;

        int result = MPI_Comm_split(temp_comm, color, 0, &comm);
        if (result != MPI_SUCCESS) throw std::runtime_error(
                "MPI error in MPI_Comm_split");
        if (ranks.size() > 0) {
            error = MPI_Comm_free(&temp_comm);
            if (error != MPI_SUCCESS) {
                throw std::runtime_error(
                        "MPI error in Commxx(MPI_Comm_free(temp_comm))");
            }
        }
    } else {
        comm = temp_comm;
    }
}
#endif

const Commxx Commxx::World(comm_type::world);
const Commxx Commxx::Null(comm_type::null);

Commxx::Commxx(comm_type type) 
    : comm(type==comm_type::null ? nullptr : new MPI_Comm(MPI_COMM_WORLD))
    , parent_comm()
    , type(type)
    , color(0)
    , key(0)
{
    if (type == comm_type::regular)
        throw std::runtime_error("only null or mpi_comm_world can be created");
}

void Commxx::construct()
{
    if (!parent_comm || parent_comm->is_null())
        throw std::runtime_error("invalid parent communicator while in construct");

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

Commxx Commxx::parent() const
{
    if (!parent_comm)
        throw std::runtime_error("Invalid parent communicator");

    return *parent_comm;
}

Commxx Commxx::dup() const
{
    if (is_null()) throw std::runtime_error("dup from a null comm");
    return Commxx(shared_from_this(), 0, rank());
}

Commxx Commxx::split(int color) const
{
    if (is_null()) throw std::runtime_error("split from a null comm");
    return Commxx(shared_from_this(), color, rank());
}

Commxx Commxx::split(int color, int key) const
{
    if (is_null()) throw std::runtime_error("split from a null comm");
    return Commxx(shared_from_this(), color, key);
}

Commxx Commxx::divide(int subgroup_size) const
{
    if (is_null()) throw std::runtime_error("divide from a null comm");

    if (size() < subgroup_size) return dup();
    if (size() % subgroup_size != 0)
    {
        throw std::runtime_error(
                "Commxx::divide(): size must be divisible by subgroup_size" );
    }

    int color = rank() / subgroup_size;
    return split(color);
}

Commxx Commxx::group(std::vector<int> const & ranks) const
{
    if (is_null()) throw std::runtime_error("group from a null comm");

    int r = rank();
    int color = (std::find(ranks.begin(), ranks.end(), r) != ranks.end())
        ? 0 : MPI_UNDEFINED;

    return Commxx(shared_from_this(), color, r);

#if 0
    MPI_Group grp;
    MPI_Comm_group(MPI_Comm(*this), &grp);

    MPI_Group subgrp;
    MPI_Group_incl(grp, ranks.size(), &ranks[0], &subgrp);

    MPI_Comm subcomm;
    MPI_Comm_create(MPI_Comm(*this), subgrp, &subcomm);

    MPI_Group_free(&grp);
    MPI_Group_free(&subgrp);

    return Commxx(subcomm, comm_create_kind::take_ownership);
#endif
}

bool operator== (Commxx const & comm1, Commxx const & comm2)
{
    if (comm1.is_null() && comm2.is_null())
        return true;

    int result;
    MPI_Comm_compare((MPI_Comm)comm1, (MPI_Comm)comm2, &result);
    return result == MPI_IDENT;
}

#if 0
Commxxs
generate_subcomms(Commxx_sptr parent_sptr, int count)
{
    Commxxs retval(0);
    std::vector<std::vector<int > > ranks(distribute_1d(*parent_sptr, count));
    for (int index = 0; index < count; ++index) {
        if (index == 0) {
            retval.push_back(
                    Commxx_sptr(new Commxx(parent_sptr, ranks.at(index))));
        } else if ((ranks.at(index) == ranks.at(index - 1))) {
            retval.push_back(retval.at(index - 1));
        } else {
            retval.push_back(
                    Commxx_sptr(new Commxx(parent_sptr, ranks.at(index))));
        }
    }
    return retval;
}

Commxx_sptr
make_optimal_spc_comm(Commxx_sptr comm_sptr, int optimal_number, bool equally_spread)
{   
    int optimal_size=std::min(optimal_number, comm_sptr->get_size());
    std::vector<int > on_ranks(optimal_size);
    int start_rank;
     if (equally_spread){
        if ((comm_sptr->get_size() %  optimal_size) !=0)  
	  throw std::runtime_error("make_optimal_spc_comm, for equal_spread  the subsize is not a divider of comm size");
	start_rank=comm_sptr->get_rank()/optimal_size;
     }
     else{
       start_rank=0;
     }      
    for (int i=0; i<optimal_size;++i){
	    on_ranks[i]=i+start_rank*optimal_size;
    }	
    Commxx_sptr ret_comm_sptr(new Commxx(comm_sptr,on_ranks)); 
    return ret_comm_sptr; 
} 
#endif
