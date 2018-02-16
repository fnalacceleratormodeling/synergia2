#include "bunch_train.h"
#include "synergia/utils/parallel_utils.h"

void 
Bunch_train::find_parent_comm_sptr()
{
  try{
    if (bunches.size()>0) {
        // check if all bunches has the same parent communicator
        MPI_Comm comm_test=bunches[0]->get_comm().get_parent_sptr()->get();
	for (size_t i = 1; i < bunches.size(); ++i) {
 	      int result;
 	      MPI_Comm_compare( comm_test, bunches[i]->get_comm().get_parent_sptr()->get(), &result);
 	      if (result != MPI_IDENT) {
             	throw std::runtime_error("Bunch_train, find_parent_comm_sptr: bunches have different parrent comunicator");
 	      } 
 	} 
 	parent_comm_sptr=bunches[0]->get_comm().get_parent_sptr();
    }
    else{
	 throw std::runtime_error(
	    "Bunch_train, find_parent_comm_sptr: number of bunches is zero, there is no commuicator");
    }   
    has_parent_comm=true;
  }
  catch (std::exception const& e) {
         std::cout<<e.what()<<std::endl;
         MPI_Abort(MPI_COMM_WORLD, 333);
  }  
}  

void 
Bunch_train::calculates_counts_and_offsets_for_impedance()
{
   
  try{
     if (!has_parent_comm) find_parent_comm_sptr(); 
     int size_parent=parent_comm_sptr->get_size();
     proc_counts_imped.resize(size_parent);
     proc_offsets_imped.resize(size_parent); 
 	  counts_and_offsets_for_impedance(*parent_comm_sptr, bunches.size(), proc_offsets_imped, proc_counts_imped);
  }
  catch (std::exception const& e) {
        std::cout<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD, 333);
  }  
}  


Commxx_sptr 
Bunch_train::get_parent_comm_sptr()
{
  if (!has_parent_comm) find_parent_comm_sptr(); 
  return parent_comm_sptr;
}  

void
Bunch_train::set_bucket_indices()
{
    std::list<int > found_indices;
    for (size_t i = 0; i < bunches.size(); ++i) {
        if (!bunches[i]->is_bucket_index_assigned()) {
            bunches[i]->set_bucket_index(i);
        }
        for (std::list<int >::const_iterator it = found_indices.begin();
                it != found_indices.end(); ++it) {
            if (*it >= bunches[i]->get_bucket_index()) {
                throw std::runtime_error(
                        "Bunch_train: bunch bucket indices must be either in strictly increasing order or all zero; otherwise wake field does not work");
            }
        }
        found_indices.push_back(bunches[i]->get_bucket_index());
    }
}

Bunch_train::Bunch_train(Bunches const& bunches, double spacing) :
                bunches(bunches),
                spacings(std::vector<double >(bunches.size() - 1, spacing)),
                 has_parent_comm(false)
{
    set_bucket_indices();
    calculates_counts_and_offsets_for_impedance();
}

Bunch_train::Bunch_train(Bunches const& bunches,
        std::vector<double > const& spacings) :
                bunches(bunches),
                spacings(spacings),
                has_parent_comm(false)
{
    if (spacings.size() != bunches.size() - 1) {
        throw std::runtime_error(
                "Bunch_train:: spacings must have length (length(bunches)-1)");
    }
    set_bucket_indices();
    calculates_counts_and_offsets_for_impedance();
}

Bunch_train::Bunch_train()
{
}

size_t
Bunch_train::get_size() const
{
    return bunches.size();
}

Bunches &
Bunch_train::get_bunches()
{
    return bunches;
}

std::vector<double > &
Bunch_train::get_spacings()
{
    return spacings;
}


std::vector< int> &
Bunch_train::get_proc_counts_for_impedance() 
{
  return proc_counts_imped;
}

std::vector< int> &
Bunch_train::get_proc_offsets_for_impedance() 
{
  return proc_offsets_imped;
}


void
Bunch_train::update_bunch_total_num()
{
    const size_t nb = get_size();
    if (nb == 0) return;

    std::vector<int> nums(nb);

    for (int i=0; i<nb; ++i)
    {
        nums[i] = bunches[i]->get_local_num();
    }

    MPI_Allreduce(MPI_IN_PLACE, &nums[0], nb, MPI_INT, MPI_SUM, 
            get_parent_comm_sptr()->get());

    for (int i=0; i<nb; ++i)
    {
        bunches[i]->set_total_num(nums[i]);
    }
}


template<class Archive>
    void
    Bunch_train::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(bunches);
        ar & BOOST_SERIALIZATION_NVP(spacings);
        ar & BOOST_SERIALIZATION_NVP(has_parent_comm);
        ar & BOOST_SERIALIZATION_NVP(parent_comm_sptr);
        ar & BOOST_SERIALIZATION_NVP(proc_counts_imped);
        ar & BOOST_SERIALIZATION_NVP(proc_offsets_imped);
    }

template
void
Bunch_train::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Bunch_train::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

Bunch_train::~Bunch_train()
{
}
