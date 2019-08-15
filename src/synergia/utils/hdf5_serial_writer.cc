#include "hdf5_serial_writer.h"



void Hdf5_serial_writer::do_setup(std::vector<int> const& data_dims)
{
    std::vector<hsize_t> chunk_dims(data_rank + 1);

    dims.resize(data_rank + 1);
    max_dims.resize(data_rank + 1);
    size.resize(data_rank + 1);
    offset.resize(data_rank + 1);

    for (int i = 0; i < data_rank; ++i) 
    {
        dims[i] = data_dims.at(i);
        max_dims[i] = data_dims.at(i);
        size[i] = data_dims.at(i);
        chunk_dims[i] = data_dims.at(i);
        offset[i] = 0;
    }

    max_dims[data_rank] = H5S_UNLIMITED;

    if (resume) 
    {
        dataset = H5Dopen(file_ptr, name.c_str(), H5P_DEFAULT);
        Hdf5_handler dataspace = H5Dget_space(dataset);

        if ( H5Sget_simple_extent_ndims(dataspace) != data_rank + 1 ) 
        {
            throw std::runtime_error(
                    "Hdf5_serial_writer::resumed data has wrong rank");
        }

        herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
        if (res < 0) throw Hdf5_exception();

        size[data_rank] = dims[data_rank];
        offset[data_rank] = dims[data_rank];
        dims[data_rank] = 1;
    } 
    else 
    {
        size[data_rank] = 0;
        offset[data_rank] = 0;
        dims[data_rank] = 1;

        if (data_size == 0) 
        {
            throw std::runtime_error(
                    "Hdf5_serial_writer: zero data size encountered");
        }

        const size_t good_chunk_size = 8192; // pulled out of air

        chunk_dims[data_rank] = 
            (data_size < good_chunk_size) ? good_chunk_size/data_size : 1;

        Hdf5_handler cparms = H5Pcreate(H5P_DATASET_CREATE);
        herr_t res = H5Pset_chunk(cparms, data_rank+1, &chunk_dims[0]);
        if (res < 0) throw Hdf5_exception();

        Hdf5_handler dataspace = H5Screate_simple(data_rank+1, &dims[0], &max_dims[0]);
        dataset = H5Dcreate(file_ptr, name.c_str(), atomic_type, 
                dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
    }

    have_setup = true;
}

void Hdf5_serial_writer::do_append(void* ptr)
{
    Hdf5_handler dataspace = H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);
    ++size[data_rank];

    herr_t res = H5Dextend(dataset, &size[0]);
    if (res < 0) throw Hdf5_exception();

    Hdf5_handler filespace = H5Dget_space(dataset);

    res = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception();

    res = H5Dwrite(dataset, atomic_type, dataspace, filespace, H5P_DEFAULT, ptr);
    if (res < 0) throw Hdf5_exception();

    ++offset[data_rank];
}


#if 0
template<>
void Hdf5_serial_writer::append<karray1d>(karray1d const& data)
{
    if (!have_setup) 
    {
        data_rank = 1;
        atomic_type = hdf5_atomic_data_type<double>();
        data_size = sizeof(double) * data.num_elements();

        std::vector<int> data_dims(data_rank);
        for (int i=0; i<data_rank; ++i) data_dims[i] = data.shape()[i];

        setup(data_dims);
    }

    Hdf5_handler dataspace = H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);
    ++size[data_rank];

    herr_t res = H5Dextend(dataset, &size[0]);
    if (res < 0) throw Hdf5_exception();

    Hdf5_handler filespace = H5Dget_space(dataset);

    res = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception();

    res = H5Dwrite(dataset, atomic_type, dataspace, filespace, H5P_DEFAULT, data.origin());
    if (res < 0) throw Hdf5_exception();

    ++offset[data_rank];
}

template<>
    void
    Hdf5_serial_writer<MArray2d_ref >::append(MArray2d_ref & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(data_rank);
            for (int i = 0; i < data_rank; ++i) {
                data_dims.at(i) = data.shape()[i];
            }
            data_size = sizeof(double) * data.num_elements();
            setup(data_dims);
        }

        Hdf5_handler dataspace = H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];

        herr_t res = H5Dextend(dataset, &size[0]);
        if (res < 0) throw Hdf5_exception();

        Hdf5_handler filespace = H5Dget_space(dataset);

        res = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &dims[0], NULL);
        if (res < 0) throw Hdf5_exception();

        res = H5Dwrite(dataset, atomic_type, dataspace, filespace, H5P_DEFAULT, data.origin());
        if (res < 0) throw Hdf5_exception();

        ++offset[data_rank];
    }

template<>
    void
    Hdf5_serial_writer<MArray3d_ref >::append(MArray3d_ref & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(data_rank);
            for (int i = 0; i < data_rank; ++i) {
                data_dims.at(i) = data.shape()[i];
            }
            data_size = sizeof(double) * data.num_elements();
            setup(data_dims);
        }

        Hdf5_handler dataspace = H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];

        herr_t res = H5Dextend(dataset, &size[0]);
        if (res < 0) throw Hdf5_exception();

        Hdf5_handler filespace = H5Dget_space(dataset);

        res = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &dims[0], NULL);
        if (res < 0) throw Hdf5_exception();

        res = H5Dwrite(dataset, atomic_type, dataspace, filespace, H5P_DEFAULT, data.origin());
        if (res < 0) throw Hdf5_exception();

        ++offset[data_rank];
    }
#endif
