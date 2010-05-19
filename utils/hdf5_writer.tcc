#include "utils/multi_array_typedefs.h"
#include <stdexcept>
#include <iostream>
inline void
h5_error_check(hid_t status)
{
    if (status != 0) {
        throw(std::runtime_error("hdf5 error"));
    }
}

// h5_atomic_typename is a local function. The generic (T) version of the
// template is undefined; only versions with specializations will compile.
template<typename T>
    inline hid_t
    h5_atomic_typename();

template<>
    inline hid_t
    h5_atomic_typename<int > ()
    {
        return H5T_NATIVE_INT;
    }

template<>
    inline hid_t
    h5_atomic_typename<double > ()
    {
        return H5T_NATIVE_DOUBLE;
    }

template<typename T>
    void
    Hdf5_writer<T >::setup(std::vector<int > const& data_dims,
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
        h5_error_check(status);
        dataset = H5Dcreate(file, name.c_str(), h5_atomic_type, dataspace,
                H5P_DEFAULT, cparms, H5P_DEFAULT);
        have_setup = true;
    }

// this is the generic case, it will only work for atomic T's that have
// h5_atomic_typename<T>() defined
template<typename T>
    Hdf5_writer<T >::Hdf5_writer(hid_t & file, std::string const& name) :
        file(file), name(name), data_rank(0), have_setup(false)
    {
    }

template<>
    Hdf5_writer<MArray1d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name);
template<>
    Hdf5_writer<MArray2d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name);

template<>
    Hdf5_writer<MArray3d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name);

template<typename T>
    void
    Hdf5_writer<T >::append(T & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(1); // dummy variable -- length really should
            // be 0, but that would not compile
            setup(data_dims, h5_atomic_typename<T > ());
        }
        ++size[data_rank];
        status = H5Dextend(dataset, &size[0]);
        h5_error_check(status);

        filespace = H5Dget_space(dataset);
        have_filespace = true;
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0],
                NULL, &dims[0], NULL);
        h5_error_check(status);
        status = H5Dwrite(dataset, h5_atomic_type, dataspace, filespace,
                H5P_DEFAULT, &data);
        h5_error_check(status);
        ++offset[data_rank];
    }

template<>
    void
    Hdf5_writer<MArray1d_ref >::append(MArray1d_ref & data);

template<>
    void
    Hdf5_writer<MArray2d_ref >::append(MArray2d_ref & data);

template<>
    void
    Hdf5_writer<MArray3d_ref >::append(MArray3d_ref & data);

template<typename T>
    Hdf5_writer<T >::~Hdf5_writer()
    {
        if (have_setup) {
            status = H5Pclose(cparms);
            h5_error_check(status);
            status = H5Dclose(dataset);
            h5_error_check(status);
            status = H5Sclose(dataspace);
            h5_error_check(status);
            if (have_filespace) {
                status = H5Sclose(filespace);
                h5_error_check(status);
            }
        }
    }
