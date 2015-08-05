#include "ff_element.h"

int FF_element::default_steps = 6; // default yoshida steps
int FF_element::default_order = 4; // default yoshida order

template<class Archive>
    void FF_element::serialize(Archive & ar, const unsigned int version)
    {
    }

template
void FF_element::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void FF_element::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void FF_element::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void FF_element::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);

FF_element::~FF_element()
{

}

