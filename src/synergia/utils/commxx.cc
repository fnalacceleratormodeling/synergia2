#include "commxx.h"
#include "parallel_utils.h"
#include <stdexcept>
#include <climits>

// hash is a local function
static size_t
hash(const char * s)
{
    size_t h = 37062913;
    while (*s)
        h = h * 101 + (unsigned char) *s++;
    return h;
}

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

Commxx::Commxx() :
        comm(MPI_COMM_WORLD ), per_host(false), ranks(0), parent_sptr(), has_this_rank_(
                true)
{
}

Commxx::Commxx(bool per_host) :
        per_host(per_host), ranks(0), parent_sptr()
{
    construct(MPI_COMM_WORLD );
}

Commxx::Commxx(Commxx_sptr parent_sptr, bool per_host) :
        per_host(per_host), ranks(0), parent_sptr()
{
    construct(parent_sptr->get());
}

Commxx::Commxx(Commxx_sptr parent_sptr, std::vector<int > const& ranks,
        bool per_host) :
        per_host(per_host), ranks(ranks), parent_sptr(parent_sptr)
{
    construct(parent_sptr->get());
}

Commxx_sptr
Commxx::get_parent_sptr() const
{
  return parent_sptr;
}  
int
Commxx::get_rank() const
{
    int error, rank;
    error = MPI_Comm_rank(comm, &rank);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_rank");
    }
    return rank;
}

int
Commxx::get_size() const
{
    int error, size;
    error = MPI_Comm_size(comm, &size);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_size");
    }
    return size;
}

bool
Commxx::has_this_rank() const
{
    return has_this_rank_;
}

MPI_Comm
Commxx::get() const
{
    return comm;
}

template<class Archive>
    void
    Commxx::save(Archive & ar, const unsigned int version) const
    {
        ar & BOOST_SERIALIZATION_NVP(per_host);
        ar & BOOST_SERIALIZATION_NVP(ranks);
        ar & BOOST_SERIALIZATION_NVP(parent_sptr);
        ar & BOOST_SERIALIZATION_NVP(has_this_rank_);
    }

template<class Archive>
    void
    Commxx::load(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(per_host);
        ar & BOOST_SERIALIZATION_NVP(ranks);
        ar & BOOST_SERIALIZATION_NVP(parent_sptr);
        ar & BOOST_SERIALIZATION_NVP(has_this_rank_);
        if (parent_sptr) {
            construct(parent_sptr->get());
        } else {
            construct(MPI_COMM_WORLD);
        }
    }

template
void
Commxx::save<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version) const;
template
void
Commxx::save<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version) const;

template
void
Commxx::load<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
template
void
Commxx::load<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Commxx::~Commxx()
{
    if (((ranks.size() > 0) || per_host) && has_this_rank_) {
        int error = MPI_Comm_free(&comm);
        if (error != MPI_SUCCESS) {
            throw std::runtime_error("MPI error in Commxx(MPI_Comm_free)");
        }
    }
}

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
