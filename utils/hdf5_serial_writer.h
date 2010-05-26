#ifndef HDF5_SERIAL_WRITER_H_
#define HDF5_SERIAL_WRITER_H_
#include <vector>
#include <string>
#include "hdf5.h"

void
h5_error_check(hid_t status);

template<typename T>
    class Hdf5_serial_writer
    {
    private:
        std::vector<hsize_t > dims, max_dims, size, offset, chunk_dims;
        herr_t status;
        int data_rank;
        std::string name;
        hid_t file, dataspace, cparms, dataset, filespace;
        hid_t h5_atomic_type;
        bool have_filespace;
        bool have_setup;
        void
        setup(std::vector<int > const& data_dims, hid_t const& file);
    public:
        Hdf5_serial_writer(hid_t & file, std::string const& name);
        void
        append(T & data);
        ~Hdf5_serial_writer();
    };

#include "utils/hdf5_serial_writer.tcc"
#endif /* HDF5_SERIAL_WRITER_H_ */
