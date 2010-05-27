#include "hdf5_chunked_array2d_writer.h"
#include "hdf5_misc.h"
#include <stdexcept>
#include <iostream>

Hdf5_chunked_array2d_writer::Hdf5_chunked_array2d_writer(hid_t & file,
        std::string const& name, Const_MArray2d_ref const & initial_data) :
    dims(2), max_dims(2), size(2), offset(2), chunk_dims(2)
{
    for (int i = 0; i < 2; ++i) {
        dims[i] = initial_data.shape()[i];
        max_dims[i] = initial_data.shape()[i];
        chunk_dims[i] = initial_data.shape()[i];
        offset[i] = 0;
    }
    size[0] = initial_data.shape()[0];
    size[1] = 0;
    max_dims[1] = H5S_UNLIMITED;
    dataspace = H5Screate_simple(2, &dims[0], &max_dims[0]);
    cparms = H5Pcreate(H5P_DATASET_CREATE);
    herr_t status = H5Pset_chunk(cparms, 2, &chunk_dims[0]);
    hdf5_error_check(status);
    dataset = H5Dcreate(file, name.c_str(), hdf5_atomic_typename<double > (),
            dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    closed = false;
    have_filespace = false;
}

void
Hdf5_chunked_array2d_writer::write_chunk(Const_MArray2d_ref const & data)
{
    if (closed) {
        throw std::runtime_error(
                "Hdf5_chunked_array2d_writer: attempt to write_chunk to a closed writer");
    }
    chunk_dims[0] = data.shape()[0];
    chunk_dims[1] = data.shape()[1];
    herr_t status = H5Pset_chunk(cparms, 2, &chunk_dims[0]);
    hdf5_error_check(status);
    size[1] += data.shape()[1];
    std::cout << "jfa: size = " << size[0] << ", " << size[1] << std::endl;
    status = H5Dextend(dataset, &size[0]);
    hdf5_error_check(status);
    filespace = H5Dget_space(dataset);
    have_filespace = true;
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL,
            &dims[0], NULL);
    hdf5_error_check(status);
    status = H5Dwrite(dataset, hdf5_atomic_typename<double > (), dataspace,
            filespace, H5P_DEFAULT, data.origin());
    hdf5_error_check(status);
    offset[1] += data.shape()[1];
}

void
Hdf5_chunked_array2d_writer::write_chunk(MArray2d_view & data)
{
    if (closed) {
        throw std::runtime_error(
                "Hdf5_chunked_array2d_writer: attempt to write_chunk to a closed writer");
    }
    chunk_dims[0] = data.shape()[0];
    chunk_dims[1] = data.shape()[1];
    herr_t status = H5Pset_chunk(cparms, 2, &chunk_dims[0]);
    hdf5_error_check(status);
    size[1] += data.shape()[1];
    std::cout << "jfa: size = " << size[0] << ", " << size[1] << std::endl;
    status = H5Dextend(dataset, &size[0]);
    hdf5_error_check(status);
    filespace = H5Dget_space(dataset);
    have_filespace = true;
    status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL,
            &dims[0], NULL);
    hdf5_error_check(status);
    status = H5Dwrite(dataset, hdf5_atomic_typename<double > (), dataspace,
            filespace, H5P_DEFAULT, data.origin());
    hdf5_error_check(status);
    offset[1] += data.shape()[1];
}

void
Hdf5_chunked_array2d_writer::close()
{
    if (!closed) {
        herr_t status = H5Pclose(cparms);
        hdf5_error_check(status);
        status = H5Dclose(dataset);
        hdf5_error_check(status);
        status = H5Sclose(dataspace);
        hdf5_error_check(status);
        if (have_filespace) {
            status = H5Sclose(filespace);
            hdf5_error_check(status);
        }
    }
    closed = true;
}

Hdf5_chunked_array2d_writer::~Hdf5_chunked_array2d_writer()
{
    close();
}
