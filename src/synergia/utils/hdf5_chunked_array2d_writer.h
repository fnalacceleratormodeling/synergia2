#ifndef HDF5_CHUNKED_ARRAY2D_WRITER_H_
#define HDF5_CHUNKED_ARRAY2D_WRITER_H_
#include <vector>
#include <string>

#include "hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"

class Hdf5_chunked_array2d_writer
{
private:
    std::vector<hsize_t > dims, max_dims, size, offset, chunk_dims;
    std::string name;

    Hdf5_handler dataset;
    Hdf5_handler atomic_type;

public:
    Hdf5_chunked_array2d_writer(hid_t file_ptr, std::string const& name,
            Const_MArray2d_view const & initial_data, int chunk_size = 0);
    Hdf5_chunked_array2d_writer(hid_t file_ptr, std::string const& name,
            Const_MArray2d_ref const & initial_data, int chunk_size = 0);
    void
    write_chunk(Const_MArray2d_ref const & data);
    void
    write_chunk(Const_MArray2d_view const & data);
    ~Hdf5_chunked_array2d_writer();
};

#endif /* HDF5_CHUNKED_ARRAY2D_WRITER_H_ */