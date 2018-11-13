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
