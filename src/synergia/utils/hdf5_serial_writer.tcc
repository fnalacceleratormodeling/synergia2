#include "utils/multi_array_typedefs.h"
#include "utils/hdf5_misc.h"
#include <stdexcept>
#include <iostream>

template<typename T>
    void
    Hdf5_serial_writer<T >::setup(std::vector<int > const& data_dims,
            hid_t const& h5_atomic_type)
    {
        this->h5_atomic_type = h5_atomic_type;
        dims.resize(data_rank + 1);
        max_dims.resize(data_rank + 1);
        size.resize(data_rank + 1);
        offset.resize(data_rank + 1);
        chunk_dims.resize(data_rank + 1);
        have_filespace = false;
        for (int i = 0; i < data_rank; ++i) {
            dims[i] = data_dims.at(i);
            max_dims[i] = data_dims.at(i);
            size[i] = data_dims.at(i);
            chunk_dims[i] = data_dims.at(i);
            offset[i] = 0;
        }
        dims[data_rank] = 1;
        max_dims[data_rank] = H5S_UNLIMITED;
        size[data_rank] = 0;
        offset[data_rank] = 0;
        chunk_dims[data_rank] = 1;
        dataspace = H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);
        cparms = H5Pcreate(H5P_DATASET_CREATE);
        status = H5Pset_chunk(cparms, data_rank + 1, &chunk_dims[0]);
        hdf5_error_check(status);
        dataset = H5Dcreate(file, name.c_str(), h5_atomic_type, dataspace,
                H5P_DEFAULT, cparms, H5P_DEFAULT);
        have_setup = true;
    }

// this is the generic case, it will only work for atomic T's that have
// hdf5_atomic_typename<T>() defined
template<typename T>
    Hdf5_serial_writer<T >::Hdf5_serial_writer(hid_t & file,
            std::string const& name) :
        data_rank(0), name(name), file(file), have_setup(false)
    {
    }

template<>
    Hdf5_serial_writer<MArray1d_ref >::Hdf5_serial_writer(hid_t & file,
            std::string const& name);
template<>
    Hdf5_serial_writer<MArray2d_ref >::Hdf5_serial_writer(hid_t & file,
            std::string const& name);

template<>
    Hdf5_serial_writer<MArray3d_ref >::Hdf5_serial_writer(hid_t & file,
            std::string const& name);

template<typename T>
    void
    Hdf5_serial_writer<T >::append(T & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(1); // dummy variable -- length really should
            // be 0, but that would not compile
            setup(data_dims, hdf5_atomic_typename<T > ());
        }
        ++size[data_rank];
        status = H5Dextend(dataset, &size[0]);
        hdf5_error_check(status);

        filespace = H5Dget_space(dataset);
        have_filespace = true;
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0],
                NULL, &dims[0], NULL);
        hdf5_error_check(status);
        status = H5Dwrite(dataset, h5_atomic_type, dataspace, filespace,
                H5P_DEFAULT, &data);
        hdf5_error_check(status);
        ++offset[data_rank];
    }

template<>
    void
    Hdf5_serial_writer<MArray1d_ref >::append(MArray1d_ref & data);

template<>
    void
    Hdf5_serial_writer<MArray2d_ref >::append(MArray2d_ref & data);

template<>
    void
    Hdf5_serial_writer<MArray3d_ref >::append(MArray3d_ref & data);

template<typename T>
    Hdf5_serial_writer<T >::~Hdf5_serial_writer()
    {
        if (have_setup) {
            status = H5Pclose(cparms);
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
    }
