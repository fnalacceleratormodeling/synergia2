#ifndef HDF5_WRITER_H_
#define HDF5_WRITER_H_
#include <vector>
#include <string>

//#include "H5Cpp.h"
#include "hdf5.h"

template<typename T>
    class Hdf5_writer
    {
    private:
        int data_rank;
        std::vector<hsize_t > dims;
        std::string name;
#if 0
        H5::H5File * file_ptr;
        H5::DataType atomic_type;
#endif
        hid_t file_ptr;
        hid_t atomic_type;

        void
        update_dims(T const& data);
        const void *
        get_data_ptr(T const& data);
    public:
        //Hdf5_writer(H5::H5File * file_ptr, std::string const& name);
        Hdf5_writer(hid_t file_ptr, std::string const& name);
        void
        write(T const & data);
        ~Hdf5_writer();
    };

#include "synergia/utils/hdf5_writer.tcc"
#endif /* HDF5_WRITER_H_ */
