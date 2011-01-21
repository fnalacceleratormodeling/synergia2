#ifndef HDF5_CHUNKED_ARRAY2D_WRITER_H_
#define HDF5_CHUNKED_ARRAY2D_WRITER_H_
#include <vector>
#include <string>
#include "hdf5.h"
#include "synergia/utils/multi_array_typedefs.h"

class Hdf5_chunked_array2d_writer
{
private:
    std::vector<hsize_t > dims, max_dims, size, offset, chunk_dims;
    std::string name;
    hid_t file, dataspace, cparms, dataset, filespace;
    bool closed, have_filespace;
public:
    Hdf5_chunked_array2d_writer(hid_t & file, std::string const& name,
            Const_MArray2d_view const & initial_data);
    Hdf5_chunked_array2d_writer(hid_t & file, std::string const& name,
            Const_MArray2d_ref const & initial_data);
    void
    write_chunk(Const_MArray2d_ref const & data);
    void
    write_chunk(Const_MArray2d_view const & data);
    void
    close();
    ~Hdf5_chunked_array2d_writer();
};

#endif /* HDF5_CHUNKED_ARRAY2D_WRITER_H_ */
