#include "commxx.h"
#include "parallel_utils.h"
#include <stdexcept>
#include <climits>


const Commxx Commxx::world;


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


Commxx::Commxx(std::shared_ptr<MPI_Comm> c)
    : comm(c)
{
}

Commxx::Commxx() 
    : comm(new MPI_Comm(MPI_COMM_WORLD))
{
}

Commxx::Commxx(MPI_Comm const & mpi_comm, comm_create_kind kind)
{
    if (mpi_comm == MPI_COMM_NULL) return;

    switch(kind)
    {
        case comm_create_kind::attach:
            comm.reset(new MPI_Comm(mpi_comm));
            break;

        case comm_create_kind::duplicate:
            MPI_Comm newcomm;
            MPI_Comm_dup(mpi_comm, &newcomm);
            comm.reset(new MPI_Comm(newcomm), comm_free());
            break;

        case comm_create_kind::take_ownership:
            comm.reset(new MPI_Comm(mpi_comm), comm_free());
            break;
    }
}

int
Commxx::get_rank() const
{
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
    int error, size;
    error = MPI_Comm_size(*comm, &size);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_size");
    }
    return size;
}

Commxx Commxx::parent() const
{
    auto pc = parent_comm.lock();
    return Commxx(pc);
}

Commxx Commxx::dup() const
{
    return Commxx(MPI_Comm(*this), comm_create_kind::duplicate);
}

Commxx Commxx::split(int color) const
{
    return split(color, rank());
}

Commxx Commxx::split(int color, int key) const
{
    MPI_Comm newcomm;
    MPI_Comm_split(MPI_Comm(*this), color, key, &newcomm);
    return Commxx(newcomm, comm_create_kind::take_ownership);
}

Commxx Commxx::group(std::vector<int> const & ranks) const
{
    MPI_Group grp;
    MPI_Comm_group(MPI_Comm(*this), &grp);

    MPI_Group subgrp;
    MPI_Group_incl(grp, ranks.size(), &ranks[0], &subgrp);

    MPI_Comm subcomm;
    MPI_Comm_create(MPI_Comm(*this), subgrp, &subcomm);

    MPI_Group_free(&grp);
    MPI_Group_free(&subgrp);

    return Commxx(subcomm, comm_create_kind::take_ownership);
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
