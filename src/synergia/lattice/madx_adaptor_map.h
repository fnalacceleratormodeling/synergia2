#ifndef MADX_ADAPTOR_MAP_H_
#define MADX_ADAPTOR_MAP_H_

#include "synergia/lattice/element_adaptor_map.h"

class MadX_adaptor_map : public Element_adaptor_map
{
public:
    MadX_adaptor_map();
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
};
BOOST_CLASS_EXPORT_KEY(MadX_adaptor_map)
typedef boost::shared_ptr<MadX_adaptor_map > MadX_adaptor_map_sptr;

#endif /* MADX_ADAPTOR_MAP_H_ */
