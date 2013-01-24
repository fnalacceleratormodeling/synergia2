#ifndef MAD8_ADAPTOR_MAP_H_
#define MAD8_ADAPTOR_MAP_H_

#include "synergia/lattice/element_adaptor_map.h"

class Mad8_adaptor_map : public Element_adaptor_map
{
public:
    Mad8_adaptor_map();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};
BOOST_CLASS_EXPORT_KEY(Mad8_adaptor_map)
typedef boost::shared_ptr<Mad8_adaptor_map > Mad8_adaptor_map_sptr;

#endif /* MAD8_ADAPTOR_MAP_H_ */
