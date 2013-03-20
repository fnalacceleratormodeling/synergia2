#ifndef COMMXX_DIVIDER_H_
#define COMMXX_DIVIDER_H_

#include "synergia/utils/commxx.h"
#include <map>

class Commxx_divider
{
private:
    std::map<Commxx_sptr, Commxx_sptr > cache;
    int subsize;
    bool per_host;
public:
    Commxx_divider();
    Commxx_divider(int subsize, bool per_host);
    virtual Commxx_sptr
    get_commxx_sptr(Commxx_sptr const& parent);
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Commxx_divider();
};

typedef boost::shared_ptr<Commxx_divider > Commxx_divider_sptr;

#endif /* COMMXX_DIVIDER_H_ */
