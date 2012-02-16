#ifndef HDF5_SERIAL_WRITER_H_
#define HDF5_SERIAL_WRITER_H_
#include <vector>
#include <string>
#include "H5Cpp.h"

#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/serialization.h"

template<typename T>
    class Hdf5_serial_writer
    {
    private:
        std::vector<hsize_t > dims, max_dims, size, offset;
        int data_rank;
        std::string name;
        Hdf5_file_sptr file_sptr;
        H5::DataSet dataset;
        H5::DataType atomic_type;
        bool have_setup;
        bool resume;
        void
        setup(std::vector<int > const& data_dims, H5::DataType atomic_type);
    public:
        Hdf5_serial_writer(Hdf5_file_sptr file_sptr, std::string const& name,
                bool resume = false);
        // Default constructor for serialization use only
        Hdf5_serial_writer();
        void
        append(T & data);
        template<class Archive>
            void
            save(Archive & ar, const unsigned int version) const
            {
                ar & BOOST_SERIALIZATION_NVP(name);
                ar & BOOST_SERIALIZATION_NVP(offset);
                ar & BOOST_SERIALIZATION_NVP(file_sptr);
            }
        template<class Archive>
            void
            load(Archive & ar, const unsigned int version)
            {
                ar & BOOST_SERIALIZATION_NVP(name);
                ar & BOOST_SERIALIZATION_NVP(offset);
                ar & BOOST_SERIALIZATION_NVP(file_sptr);
                resume = true;
                have_setup = false;
            }
        BOOST_SERIALIZATION_SPLIT_MEMBER()
        ~Hdf5_serial_writer();
    };

#include "synergia/utils/hdf5_serial_writer.tcc"
#endif /* HDF5_SERIAL_WRITER_H_ */
