#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"
#include <stdexcept>
#include <iostream>

template<typename T>
    void
    Hdf5_serial_writer<T >::setup(std::vector<int > const& data_dims)
    {
        std::vector<hsize_t > chunk_dims(data_rank + 1);

        dims.resize(data_rank + 1);
        max_dims.resize(data_rank + 1);
        size.resize(data_rank + 1);
        offset.resize(data_rank + 1);

        for (int i = 0; i < data_rank; ++i) {
            dims[i+1] = data_dims.at(i);
            max_dims[i+1] = data_dims.at(i);
            size[i+1] = data_dims.at(i);
            chunk_dims[i+1] = data_dims.at(i);
            offset[i+1] = 0;
        }

        max_dims[0] = H5S_UNLIMITED;

        if (resume) {
            dataset = H5Dopen(file_sptr->get_h5file(), name.c_str(), H5P_DEFAULT);
            Hdf5_handler dataspace = H5Dget_space(dataset);
            int file_rank = H5Sget_simple_extent_ndims(dataspace);

            if (file_rank != data_rank + 1) {
                throw std::runtime_error(
                        "Hdf5_serial_writer::resumed data has wrong rank");
            }

            herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
            if (res < 0) throw Hdf5_exception();

            size[0] = dims[0];
            offset[0] = dims[0];
            dims[0] = 1;
        } else {
            size[0] = 0;
            offset[0] = 0;
            dims[0] = 1;

            if (data_size == 0) {
                throw std::runtime_error("Hdf5_serial_writer: zero data size encountered");
            }
            const size_t good_chunk_size = 8192; // pulled out of air
            if (data_size < good_chunk_size) {
                chunk_dims[0] = good_chunk_size/data_size;
            } else {
                chunk_dims[0] = 1;
            }

            Hdf5_handler cparms = H5Pcreate(H5P_DATASET_CREATE);
            herr_t res = H5Pset_chunk(cparms, data_rank + 1, &chunk_dims[0]);
            if (res < 0) throw Hdf5_exception();

            Hdf5_handler dataspace = H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);
            dataset = H5Dcreate(file_sptr->get_h5file(), name.c_str(), atomic_type, 
                    dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
        }

        have_setup = true;
    }

// this is the generic case, it will only work for atomic T's that have
// hdf5_atomic_typename<T>() defined
template<typename T>
    Hdf5_serial_writer<T >::Hdf5_serial_writer(Hdf5_file_sptr file_sptr,
            std::string const& name, bool resume) :
            data_rank(0), name(name), file_sptr(file_sptr), atomic_type(
                    hdf5_atomic_data_type<T >()), have_setup(false), resume(
                    resume), data_size(sizeof(T))
    {
    }

template<typename T>
    Hdf5_serial_writer<T >::Hdf5_serial_writer() : atomic_type(hdf5_atomic_data_type<T > ())
    {
    }

template<>
    Hdf5_serial_writer<MArray1d_ref >::Hdf5_serial_writer(
            Hdf5_file_sptr file_sptr, std::string const& name, bool resume);

template<>
    Hdf5_serial_writer<MArray1d_ref >::Hdf5_serial_writer();

template<>
    Hdf5_serial_writer<MArray2d_ref >::Hdf5_serial_writer(
            Hdf5_file_sptr file_sptr, std::string const& name, bool resume);

template<>
    Hdf5_serial_writer<MArray2d_ref >::Hdf5_serial_writer();

template<>
    Hdf5_serial_writer<MArray3d_ref >::Hdf5_serial_writer(
            Hdf5_file_sptr file_sptr, std::string const& name, bool resume);

template<>
    Hdf5_serial_writer<MArray3d_ref >::Hdf5_serial_writer();

template<typename T>
    void
    Hdf5_serial_writer<T >::append(T & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(1); // dummy variable -- length really should
            // be 0, but that would not compile
            setup(data_dims);
        }

        Hdf5_handler dataspace = H5Screate_simple(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[0];

        herr_t res = H5Dextend(dataset, &size[0]);
        if (res < 0) throw Hdf5_exception();

        Hdf5_handler filespace = H5Dget_space(dataset);

        res = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0], NULL, &dims[0], NULL);
        if (res < 0) throw Hdf5_exception();

        res = H5Dwrite(dataset, atomic_type, dataspace, filespace, H5P_DEFAULT, &data);
        if (res < 0) throw Hdf5_exception();

        ++offset[0];
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
    }
