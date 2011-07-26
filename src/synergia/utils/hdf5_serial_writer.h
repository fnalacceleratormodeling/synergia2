#ifndef HDF5_SERIAL_WRITER_H_
#define HDF5_SERIAL_WRITER_H_
#include <vector>
#include <string>
#include "H5Cpp.h"

template<typename T>
    class Hdf5_serial_writer
    {
    private:
        std::vector<hsize_t > dims, max_dims, size, offset, chunk_dims;
        herr_t status;
        int data_rank;
        std::string name;
        H5::H5File file;
        H5::DataSet dataset;
        H5::DataType h5_atomic_type;
        bool have_filespace;
        bool have_setup;
        void
        setup(std::vector<int > const& data_dims, H5::DataType h5_atomic_type);
    public:
        Hdf5_serial_writer(H5::H5File & file, std::string const& name);
        void
        append(T & data);
        ~Hdf5_serial_writer();
    };

#include "synergia/utils/hdf5_serial_writer.tcc"
#endif /* HDF5_SERIAL_WRITER_H_ */
