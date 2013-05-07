#ifndef COMMXX_H_
#define COMMXX_H_

#include <vector>
#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include "mpi.h"

#include "synergia/utils/serialization.h"

class Commxx;
typedef boost::shared_ptr<Commxx > Commxx_sptr; // syndoc:include

/// Commxx is a wrapper around MPI communicator (MPI_Comm) objects.
///
/// Python:  (jfa: needs update) The equivalent functionality is provided by mpi4py Comm objects.
/// mpi4py Comm objects may be passed from python anywhere a Commxx object is
/// expected.
class Commxx
{
private:
    MPI_Comm comm;
    bool per_host;
    std::vector<int > ranks;
    Commxx_sptr parent_sptr;
    bool has_this_rank_;
    void
    construct(MPI_Comm const& parent_mpi_comm);
    // This dummy constructor exists only to force compilation errors in
    // code that uses an older interface to this class.
    Commxx(MPI_Comm const& comm){}
    // prevent copying
    Commxx(Commxx const& commxx){}
    Commxx&
    operator=(Commxx const& commxx){ return *this;}
public:
    /// Construct a Commxx object using MPI_COMM_WORLD
    Commxx();

    /// Construct a Commxx object, optionally creating separate communicators on each
    /// unique host for communication avoidance
    Commxx(bool per_host);

    /// Construct a Commxx object based on the parent communicator, optionally
    /// creating separate communicators on each unique host for communication avoidance
    Commxx(Commxx_sptr parent_sptr, bool per_host);

    /// Construct a Commxx object using only the specified ranks on the parent
    /// communicator
    Commxx(Commxx_sptr parent_sptr, std::vector<int > const& ranks,
            bool per_host = false);

	    
    Commxx_sptr
    get_parent_sptr() const;    
    /// Get communicator rank
    int
    get_rank() const;

    /// Get communicator size
    int
    get_size() const;

    /// Test to see if the communicator contains this rank
    bool
    has_this_rank() const;

    /// Extract the MPI_comm object wrapped by the Commxx instance.
    MPI_Comm
    get() const;

    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const;
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version);
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    ~Commxx();
};

typedef std::vector<Commxx_sptr > Commxxs;

Commxxs
generate_subcomms(Commxx_sptr parent_sptr, int count);

Commxx_sptr
make_optimal_spc_comm(Commxx_sptr parent_sptr, int optimal_number, bool equally_spread=false);


#endif /* COMMXX_H_ */
