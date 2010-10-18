#ifndef MULTI_ARRAY_SERIALIZATION_H_
#define MULTI_ARRAY_SERIALIZATION_H_

#include "boost/multi_array.hpp"
#include "boost/serialization/nvp.hpp"
#include "boost/serialization/utility.hpp"
#include "boost/serialization/version.hpp"

// n.b. Serialization of 3d arrays not yet implemented.

namespace boost
{
    namespace serialization
    {
        template<typename T, class Archive>
            void
            load(Archive & ar, boost::multi_array<T, 1 > & t,
                    const unsigned int file_version)
            {
                typedef boost::multi_array<T, 1 > multi_array_;
                typename multi_array_::size_type shape_0;
                ar >> BOOST_SERIALIZATION_NVP(shape_0);
                t.resize(boost::extents[shape_0]);
                ar >> make_array(t.data(), t.num_elements());
            }
        // jfa: this code can't work as written...
        //        template<typename T, class Archive>
        //            void
        //            load(Archive & ar, boost::multi_array_ref<T, 1 > & t,
        //                    const unsigned int file_version)
        //            {
        //                typedef boost::multi_array<T, 1 > multi_array_;
        //                typename multi_array_::size_type shape_0;
        //                ar >> BOOST_SERIALIZATION_NVP(shape_0);
        //                t.resize(boost::extents[shape_0]);
        //                ar >> make_array(t.data(), t.num_elements());
        //            }
        template<typename T, class Archive>
            void
            load(Archive & ar, boost::multi_array<T, 2 > & t,
                    const unsigned int file_version)
            {
                typedef boost::multi_array<T, 2 > multi_array_;
                typename multi_array_::size_type shape_0, shape_1;
                ar >> BOOST_SERIALIZATION_NVP(shape_0);
                ar >> BOOST_SERIALIZATION_NVP(shape_1);
                t.resize(boost::extents[shape_0][shape_1]);
                ar >> make_array(t.data(), t.num_elements());
            }
        template<typename T, typename Archive>
            void
            save(Archive & ar, const boost::multi_array<T, 1 > & t,
                    const unsigned int file_version)
            {
                typedef boost::multi_array<T, 1 > multi_array_;
                typename multi_array_::size_type shape_0 = t.shape()[0];
                ar << BOOST_SERIALIZATION_NVP(shape_0);
                ar << make_array(t.data(), t.num_elements());
            }
        template<typename T, typename Archive>
            void
            save(Archive & ar, const boost::multi_array_ref<T, 1 > & t,
                    const unsigned int file_version)
            {
                typedef boost::multi_array<T, 1 > multi_array_;
                typename multi_array_::size_type shape_0 = t.shape()[0];
                ar << BOOST_SERIALIZATION_NVP(shape_0);
                ar << make_array(t.data(), t.num_elements());
            }
        template<typename T, typename Archive>
            void
            save(Archive & ar, const boost::multi_array<T, 2 > & t,
                    const unsigned int file_version)
            {
                typedef boost::multi_array<T, 2 > multi_array_;
                typename multi_array_::size_type shape_0 = t.shape()[0];
                typename multi_array_::size_type shape_1 = t.shape()[1];
                ar << BOOST_SERIALIZATION_NVP(shape_0);
                ar << BOOST_SERIALIZATION_NVP(shape_1);
                ar << make_array(t.data(), t.num_elements());
            }
        template<typename T, typename Archive>
            void
            save(Archive & ar, const boost::multi_array_ref<T, 2 > & t,
                    const unsigned int file_version)
            {
                typedef boost::multi_array<T, 2 > multi_array_;
                typename multi_array_::size_type shape_0 = t.shape()[0];
                typename multi_array_::size_type shape_1 = t.shape()[1];
                ar << BOOST_SERIALIZATION_NVP(shape_0);
                ar << BOOST_SERIALIZATION_NVP(shape_1);
                ar << make_array(t.data(), t.num_elements());
            }
        template<typename T, class Archive>
            void
            serialize(Archive & ar, boost::multi_array<T, 1 > & t,
                    const unsigned int file_version)
            {
                split_free(ar, t, file_version);
            }
        template<typename T, class Archive>
            void
            serialize(Archive & ar, boost::multi_array_ref<T, 1 > & t,
                    const unsigned int file_version)
            {
                split_free(ar, t, file_version);
            }
        template<typename T, class Archive>
            void
            serialize(Archive & ar, boost::multi_array<T, 2 > & t,
                    const unsigned int file_version)
            {
                split_free(ar, t, file_version);
            }
        template<typename T, class Archive>
            void
            serialize(Archive & ar, boost::multi_array_ref<T, 2 > & t,
                    const unsigned int file_version)
            {
                split_free(ar, t, file_version);
            }
    }
}

#endif /* MULTI_ARRAY_SERIALIZATION_H_ */
