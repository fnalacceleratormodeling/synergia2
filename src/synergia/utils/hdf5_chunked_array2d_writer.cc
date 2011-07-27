#include "hdf5_chunked_array2d_writer.h"
#include "hdf5_misc.h"
#include <stdexcept>
#include <iostream>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

Hdf5_chunked_array2d_writer::Hdf5_chunked_array2d_writer(H5File & file,
        std::string const& name, Const_MArray2d_view const & initial_data) :
    file(file), dims(2), max_dims(2), size(2), offset(2), chunk_dims(2),
            atomic_type(hdf5_atomic_data_type<double > ())
{
    for (int i = 0; i < 2; ++i) {
        dims[i] = initial_data.shape()[i];
        max_dims[i] = initial_data.shape()[i];
        chunk_dims[i] = initial_data.shape()[i];
        offset[i] = 0;
    }
    size[0] = 0;
    size[1] = initial_data.shape()[1];
    max_dims[0] = H5S_UNLIMITED;
    DSetCreatPropList cparms;
    cparms.setChunk(2, &chunk_dims[0]);
    DataSpace dataspace(2, &dims[0], &max_dims[0]);
    dataset = file.createDataSet(name.c_str(), atomic_type, dataspace, cparms);
    closed = false;
    have_filespace = false;
}

Hdf5_chunked_array2d_writer::Hdf5_chunked_array2d_writer(H5File & file,
        std::string const& name, Const_MArray2d_ref const & initial_data) :
    file(file), dims(2), max_dims(2), size(2), offset(2), chunk_dims(2),
            atomic_type(hdf5_atomic_data_type<double > ())
{
    for (int i = 0; i < 2; ++i) {
        dims[i] = initial_data.shape()[i];
        max_dims[i] = initial_data.shape()[i];
        chunk_dims[i] = initial_data.shape()[i];
        offset[i] = 0;
    }
    size[0] = 0;
    size[1] = initial_data.shape()[1];
    max_dims[0] = H5S_UNLIMITED;
    DSetCreatPropList cparms;
    cparms.setChunk(2, &chunk_dims[0]);
    DataSpace dataspace(2, &dims[0], &max_dims[0]);
    dataset = file.createDataSet(name.c_str(), atomic_type, dataspace, cparms);
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
    DSetCreatPropList cparms;
    cparms.setChunk(2, &chunk_dims[0]);
    size[0] += data.shape()[0];
    dataset.extend(&size[0]);
    DataSpace filespace = dataset.getSpace();
    have_filespace = true;
    filespace.selectHyperslab(H5S_SELECT_SET, &chunk_dims[0], &offset[0]);
    DataSpace dataspace(2, &chunk_dims[0], &max_dims[0]);
    dataset.write(data.origin(), atomic_type, dataspace, filespace);
    offset[0] += data.shape()[0];
}

void
Hdf5_chunked_array2d_writer::write_chunk(Const_MArray2d_view const & data)
{
    if (closed) {
        throw std::runtime_error(
                "Hdf5_chunked_array2d_writer: attempt to write_chunk to a closed writer");
    }
    chunk_dims[0] = data.shape()[0];
    chunk_dims[1] = data.shape()[1];
    DSetCreatPropList cparms;
    cparms.setChunk(2, &chunk_dims[0]);
    size[0] += data.shape()[0];
    dataset.extend(&size[0]);
    DataSpace filespace = dataset.getSpace();
    have_filespace = true;
    filespace.selectHyperslab(H5S_SELECT_SET, &chunk_dims[0], &offset[0]);
    DataSpace dataspace(2, &chunk_dims[0], &max_dims[0]);
    dataset.write(data.origin(), atomic_type, dataspace, filespace);
    offset[0] += data.shape()[0];
}

void
Hdf5_chunked_array2d_writer::close()
{
    closed = true;
}

Hdf5_chunked_array2d_writer::~Hdf5_chunked_array2d_writer()
{
    close();
}
