#ifndef FF_ELEMENT_MAP_H
#define FF_ELEMENT_MAP_H

#include <map>
#include <string>

#include "synergia/libff/ff_element.h"

class FF_element_map
{
private:
    std::map<std::string, FF_element_sptr > element_map;

public:
    FF_element_map();
    bool has_adaptor(std::string const& type) const;
    void set_adaptor(std::string const& type, FF_element_sptr element_sptr);
    FF_element_sptr get_adaptor(std::string const& type) const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    ~FF_element_map();
};

#endif // FF_ELEMENT_MAP_H
