#ifndef HDF5_WRITER_H_
#define HDF5_WRITER_H_
#include <vector>
#include <string>

//#include "H5Cpp.h"
#include "hdf5.h"

#include "synergia/utils/hdf5_misc.h"

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

        // use hid_t type since the writer doesnt own the file object
        hid_t file_ptr;    
        Hdf5_handler atomic_type;

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
