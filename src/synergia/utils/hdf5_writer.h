#ifndef HDF5_WRITER_H_
#define HDF5_WRITER_H_
#include <vector>
#include <string>
#include "hdf5.h"

template<typename T>
    class Hdf5_writer
    {
    private:
        int data_rank;
        std::vector<hsize_t > dims;
        std::string name;
        hid_t file;
        hid_t atomic_type;
        void
        update_dims(T const& data);
        const void *
        get_data_ptr(T const& data);
    public:
        Hdf5_writer(hid_t & file, std::string const& name);
        void
        write(T const & data);
        ~Hdf5_writer();
    };

#include "synergia/utils/hdf5_writer.tcc"
#endif /* HDF5_WRITER_H_ */
