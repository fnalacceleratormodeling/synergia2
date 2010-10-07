#ifndef XML_SERIALIZATION_H_
#define XML_SERIALIZATION_H_

#include <fstream>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/version.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>

template<typename T>
    void
    xml_save(T const& object, const char *filename)
    {
        std::ofstream output_stream(filename);
        boost::archive::xml_oarchive output_archive(output_stream);
        output_archive << BOOST_SERIALIZATION_NVP(object);
        output_stream.close();
    }

template<typename T>
    void
    xml_load(T & object, const char *filename)
    {
        std::ifstream input_stream(filename);
        boost::archive::xml_iarchive input_archive(input_stream);
        input_archive >> BOOST_SERIALIZATION_NVP(object);
        input_stream.close();
    }

#endif /* XML_SERIALIZATION_H_ */
