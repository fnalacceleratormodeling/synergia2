#ifndef ELEMENT_ADAPTOR_MAP_H_
#define ELEMENT_ADAPTOR_MAP_H_

#include "synergia/lattice/element_adaptor.h"

class Element_adaptor_map
{
private:
    std::map<std::string, Element_adaptor_sptr > adaptor_map;

public:
    Element_adaptor_map();
    virtual void
    set_adaptor(std::string const& name,
            Element_adaptor_sptr element_adaptor_sptr);
    virtual bool
    has_adaptor(std::string const& name) const;
    virtual Element_adaptor_sptr
    get_adaptor(std::string const& name) const;
    virtual std::list<std::string >
    get_adaptor_names() const;
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version);
    virtual
    ~Element_adaptor_map() = 0;
};
BOOST_CLASS_EXPORT_KEY(Element_adaptor_map)
typedef boost::shared_ptr<Element_adaptor_map > Element_adaptor_map_sptr;

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

#endif /* ELEMENT_ADAPTOR_MAP_H_ */
