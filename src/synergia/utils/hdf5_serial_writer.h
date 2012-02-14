#ifndef HDF5_SERIAL_WRITER_H_
#define HDF5_SERIAL_WRITER_H_
#include <vector>
#include <string>
#include "H5Cpp.h"

#include "synergia/utils/hdf5_file.h"

template<typename T>
    class Hdf5_serial_writer
    {
    private:
        std::vector<hsize_t > dims, max_dims, size, offset, chunk_dims;
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
        void
        append(T & data);
        ~Hdf5_serial_writer();
    };

#include "synergia/utils/hdf5_serial_writer.tcc"
#endif /* HDF5_SERIAL_WRITER_H_ */
