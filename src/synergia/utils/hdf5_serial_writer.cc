#include "hdf5_serial_writer.h"

template<>
    Hdf5_serial_writer<MArray1d_ref >::Hdf5_serial_writer(
            Hdf5_file_sptr file_sptr, std::string const& name, bool resume) :
            data_rank(MArray1d_ref::dimensionality), name(name), file_sptr(
                    file_sptr), atomic_type(hdf5_atomic_data_type<double >()), have_setup(
                    false), resume(resume), data_size(0)
    {
    }

template<>
    Hdf5_serial_writer<MArray1d_ref >::Hdf5_serial_writer() :
            atomic_type(hdf5_atomic_data_type<double >())
    {
    }

template<>
    Hdf5_serial_writer<MArray2d_ref >::Hdf5_serial_writer(
            Hdf5_file_sptr file_sptr, std::string const& name, bool resume) :
            data_rank(MArray2d_ref::dimensionality), name(name), file_sptr(
                    file_sptr), atomic_type(hdf5_atomic_data_type<double >()), have_setup(
                    false), resume(resume), data_size(0)
    {
    }

template<>
    Hdf5_serial_writer<MArray2d_ref >::Hdf5_serial_writer() :
            atomic_type(hdf5_atomic_data_type<double >())
    {
    }

template<>
    Hdf5_serial_writer<MArray3d_ref >::Hdf5_serial_writer(
            Hdf5_file_sptr file_sptr, std::string const& name, bool resume) :
            data_rank(MArray3d_ref::dimensionality), name(name), file_sptr(
                    file_sptr), atomic_type(hdf5_atomic_data_type<double >()), have_setup(
                    false), resume(resume), data_size(0)
    {
    }

template<>
    Hdf5_serial_writer<MArray3d_ref >::Hdf5_serial_writer() :
            atomic_type(hdf5_atomic_data_type<double >())
    {
    }

template<>
    void
    Hdf5_serial_writer<MArray1d_ref >::append(MArray1d_ref & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(data_rank);
            for (int i = 0; i < data_rank; ++i) {
                data_dims.at(i) = data.shape()[i];
            }
            data_size = sizeof(double) * data.num_elements();
            setup(data_dims);
        }
        DataSpace dataspace(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];
        dataset.extend(&size[0]);

        DataSpace filespace = dataset.getSpace();
        filespace.selectHyperslab(H5S_SELECT_SET, &dims[0], &offset[0]);
        dataset.write(data.origin(), atomic_type, dataspace, filespace);
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
        DataSpace dataspace(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];
        dataset.extend(&size[0]);

        DataSpace filespace = dataset.getSpace();
        filespace.selectHyperslab(H5S_SELECT_SET, &dims[0], &offset[0]);
        dataset.write(data.origin(), atomic_type, dataspace, filespace);
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
        DataSpace dataspace(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];
        dataset.extend(&size[0]);

        DataSpace filespace = dataset.getSpace();
        filespace.selectHyperslab(H5S_SELECT_SET, &dims[0], &offset[0]);
        dataset.write(data.origin(), atomic_type, dataspace, filespace);
        ++offset[data_rank];
    }
