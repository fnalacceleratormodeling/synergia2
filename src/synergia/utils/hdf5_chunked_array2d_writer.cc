#include "hdf5_chunked_array2d_writer.h"
#include "hdf5_misc.h"
#include <stdexcept>
#include <iostream>

Hdf5_chunked_array2d_writer::Hdf5_chunked_array2d_writer(
        hid_t file_ptr,
        std::string const& name, 
        Const_MArray2d_view const & initial_data,
        int chunk_size ) 
    : dims(2), max_dims(2), size(2), offset(2), chunk_dims(2)
    , atomic_type(hdf5_atomic_data_type<double >())
{
    for (int i = 0; i < 2; ++i) 
    {
        dims[i] = initial_data.shape()[i];
        max_dims[i] = initial_data.shape()[i];
        chunk_dims[i] = initial_data.shape()[i];
        offset[i] = 0;
    }

    if (chunk_size > 0) 
    {
        chunk_dims[0] = chunk_size;
    }

    // always start with dims[0] = 0 for the initial size
    // it will be extended later when writing chunks
    dims[0] = 0;

    size[0] = 0;
    size[1] = initial_data.shape()[1];
    max_dims[0] = H5S_UNLIMITED;

    Hdf5_handler cparms = H5Pcreate(H5P_DATASET_CREATE);
    herr_t res = H5Pset_chunk(cparms, 2, &chunk_dims[0]);
    if (res < 0) throw Hdf5_exception();

    Hdf5_handler dataspace = H5Screate_simple(2, &dims[0], &max_dims[0]);
    dataset = H5Dcreate(file_ptr, name.c_str(), atomic_type,
            dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
}

Hdf5_chunked_array2d_writer::Hdf5_chunked_array2d_writer(
        hid_t file_ptr,
        std::string const& name, 
        Const_MArray2d_ref const & initial_data,
        int chunk_size ) 
    : dims(2), max_dims(2), size(2), offset(2), chunk_dims(2)
    , atomic_type(hdf5_atomic_data_type<double >())
{
    for (int i = 0; i < 2; ++i) 
    {
        dims[i] = initial_data.shape()[i];
        max_dims[i] = initial_data.shape()[i];
        chunk_dims[i] = initial_data.shape()[i];
        offset[i] = 0;
    }

    if (chunk_size > 0) 
    {
        chunk_dims[0] = chunk_size;
    }

    // always start with dims[0] = 0 for the initial size
    // it will be extended later when writing chunks
    dims[0] = 0;

    size[0] = 0;
    size[1] = initial_data.shape()[1];
    max_dims[0] = H5S_UNLIMITED;

    Hdf5_handler cparms = H5Pcreate(H5P_DATASET_CREATE);
    herr_t res = H5Pset_chunk(cparms, 2, &chunk_dims[0]);
    if (res < 0) throw Hdf5_exception();

    Hdf5_handler dataspace = H5Screate_simple(2, &dims[0], &max_dims[0]);
    dataset = H5Dcreate(file_ptr, name.c_str(), atomic_type,
            dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
}

void
Hdf5_chunked_array2d_writer::write_chunk(Const_MArray2d_ref const & data)
{
    chunk_dims[0] = data.shape()[0];
    chunk_dims[1] = data.shape()[1];

    Hdf5_handler cparms = H5Pcreate(H5P_DATASET_CREATE);
    herr_t res = H5Pset_chunk(cparms, 2, &chunk_dims[0]);
    if (res < 0) throw Hdf5_exception();

    size[0] += data.shape()[0];
    res = H5Dextend(dataset, &size[0]);
    if (res < 0) throw Hdf5_exception();

    Hdf5_handler filespace = H5Dget_space(dataset);
    res = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &chunk_dims[0], NULL);
    if (res < 0) throw Hdf5_exception();

    Hdf5_handler dataspace = H5Screate_simple(2, &chunk_dims[0], &max_dims[0]);

    // this is to make sure that the data are written in c storage order
    // since it is only writing a chunk of data from a bigger array (e.g.,
    // a chunk of particle data). cannot do that in the fortran order
    int nx = data.shape()[0];
    int ny = data.shape()[1];
    int x0 = data.index_bases()[0];
    int y0 = data.index_bases()[1];

    std::vector<double> temp(nx * ny);

    for (int i = x0; i < x0 + nx; ++i)
    {
        for (int j = y0; j < y0 + ny; ++j)
        {
            temp[i*ny + j] = data[i][j];
        }
    }

    res = H5Dwrite(dataset, atomic_type, dataspace, filespace, H5P_DEFAULT, &temp[0]);
    if (res < 0) throw Hdf5_exception();

    offset[0] += data.shape()[0];
}

void
Hdf5_chunked_array2d_writer::write_chunk(Const_MArray2d_view const & data)
{
    chunk_dims[0] = data.shape()[0];
    chunk_dims[1] = data.shape()[1];

    Hdf5_handler cparms = H5Pcreate(H5P_DATASET_CREATE);
    herr_t res = H5Pset_chunk(cparms, 2, &chunk_dims[0]);
    if (res < 0) throw Hdf5_exception();

    size[0] += data.shape()[0];
    res = H5Dextend(dataset, &size[0]);
    if (res < 0) throw Hdf5_exception();

    Hdf5_handler filespace = H5Dget_space(dataset);
    res = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &chunk_dims[0], NULL);
    if (res < 0) throw Hdf5_exception();

    Hdf5_handler dataspace = H5Screate_simple(2, &chunk_dims[0], &max_dims[0]);

    // see notes above
    int nx = data.shape()[0];
    int ny = data.shape()[1];
    int x0 = data.index_bases()[0];
    int y0 = data.index_bases()[1];

    std::vector<double> temp(nx * ny);

    for (int i = x0; i < x0 + nx; ++i)
    {
        for (int j = y0; j < y0 + ny; ++j)
        {
            temp[i*ny + j] = data[i][j];
        }
    }

    res = H5Dwrite(dataset, atomic_type, dataspace, filespace, H5P_DEFAULT, &temp[0]);
    if (res < 0) throw Hdf5_exception();

    offset[0] += data.shape()[0];
}

Hdf5_chunked_array2d_writer::~Hdf5_chunked_array2d_writer()
{
}
