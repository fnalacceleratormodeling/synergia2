#ifndef BUNCH_TRAIN_H_
#define BUNCH_TRAIN_H_

#include "synergia/bunch/bunch.h"

class Bunch_train
{
private:
    Bunches bunches;
    std::vector<double > spacings;      
    bool has_parent_comm;
    Commxx_sptr parent_comm_sptr;
    std::vector< int> proc_counts;
    std::vector< int> proc_offsets;
    void
    set_bucket_indices();
    void find_parent_comm_sptr();
    void calculates_counts_and_offsets();
public:
    Bunch_train(Bunches const& bunches, double spacing);
    Bunch_train(Bunches const& bunches, std::vector<double > const& spacings);
    // Default constructor for serialization use only
    Bunch_train();
    Commxx_sptr 
    get_parent_comm_sptr();
    size_t
    get_size() const;
    Bunches &
    get_bunches();
    std::vector<double > &
    get_spacings();
    std::vector< int> &
    get_proc_counts();
    std::vector< int> &
    get_proc_offsets();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~Bunch_train();
};

typedef boost::shared_ptr<Bunch_train > Bunch_train_sptr;

#endif /* BUNCH_TRAIN_H_ */
