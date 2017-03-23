#include "fast_mapping.h"
//#include "synergia/lattice/chef_utils.h"

//#define MANUAL_LOOP_UNROLL shockingly_yes

std::string
read_line_ignoring_comments(std::ifstream &file)
{
    std::string line;
    getline(file, line);
    if (line.size() == 0) {
        return line;
    }
    while (((line[0] == '#') || (line.size() == 0))) {
        getline(file, line);
    }
    return line;
}




template<>
template<>
void
Fast_mapping_term<double>::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);


template<>
template<>
void
Fast_mapping_term<double>::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template<>
template<>
void
Fast_mapping_term<double>::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template<>
template<>
void
Fast_mapping_term<double>::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);


template<>
template<>
void
Fast_mapping_term<std::complex<double> >::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);


template<>
template<>
void
Fast_mapping_term<std::complex<double> >::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template<>
template<>
void
Fast_mapping_term<std::complex<double> >::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template<>
template<>
void
Fast_mapping_term<std::complex<double> >::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);








template<>
template<>
void
TFast_mapping<double>::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<double>::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<double>::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<double>::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);


template<>
template<>
void
TFast_mapping<std::complex<double> >::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<std::complex<double> >::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<std::complex<double> >::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template<>
template<>
void
TFast_mapping<std::complex<double> >::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
