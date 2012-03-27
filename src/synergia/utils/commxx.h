#ifndef COMMXX_H_
#define COMMXX_H_

#include <stdexcept>
#include <boost/shared_ptr.hpp>
#include "mpi.h"

#include "synergia/utils/serialization.h"

class Commxx;
typedef boost::shared_ptr<Commxx > Commxx_sptr;

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
public:
    /// Construct a Commxx object using MPI_COMM_WORLD
    Commxx();

    /// Construct a Commxx object, optionally creating separate communicators on each
    /// unique host for communication avoidance
    Commxx(bool per_host);

    /// Construct a Commxx object using only the specified ranks on the parent
    /// communicator
    Commxx(Commxx_sptr parent_sptr, std::vector<int > const& ranks,
            bool per_host = false);

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
        save(Archive & ar, const unsigned int version) const
        {
            ar & BOOST_SERIALIZATION_NVP(per_host);
            ar & BOOST_SERIALIZATION_NVP(ranks);
            ar & BOOST_SERIALIZATION_NVP(parent_sptr);
            ar & BOOST_SERIALIZATION_NVP(has_this_rank_);
        }
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version)
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
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    ~Commxx();
};

#endif /* COMMXX_H_ */
