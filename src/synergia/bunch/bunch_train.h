#ifndef BUNCH_TRAIN_H_
#define BUNCH_TRAIN_H_

#include "synergia/bunch/bunch.h"

class Bunch_train
{
private:
    Bunches bunches;
    std::vector<double > spacings;

public:
    Bunch_train(Bunches const& bunches, double spacing);
    Bunch_train(Bunches const& bunches, std::vector<double > const& spacings);
    // Default constructor for serialization use only
    Bunch_train();
    size_t
    get_size() const;
    Bunches &
    get_bunches();
    std::vector<double > &
    get_spacings();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~Bunch_train();
};

typedef boost::shared_ptr<Bunch_train > Bunch_train_sptr;

#endif /* BUNCH_TRAIN_H_ */
