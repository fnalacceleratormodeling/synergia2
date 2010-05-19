#include "hdf5_writer.h"

template<>
    Hdf5_writer<MArray1d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name) :
        file(file), name(name), data_rank(MArray1d_ref::dimensionality),
                have_setup(false)
    {
    }

template<>
    Hdf5_writer<MArray2d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name) :
        file(file), name(name), data_rank(MArray2d_ref::dimensionality),
                have_setup(false)
    {
    }

template<>
    Hdf5_writer<MArray3d_ref >::Hdf5_writer(hid_t & file,
            std::string const& name) :
        file(file), name(name), data_rank(MArray3d_ref::dimensionality),
                have_setup(false)
    {
    }

template<>
    void
    Hdf5_writer<MArray1d_ref >::append(MArray1d_ref & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(data_rank);
            for (int i = 0; i < data_rank; ++i) {
                data_dims.at(i) = data.shape()[i];
            }
            setup(data_dims, h5_atomic_typename<double > ());
        }
        ++size[data_rank];
        status = H5Dextend(dataset, &size[0]);

        filespace = H5Dget_space(dataset);
        have_filespace = true;
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0],
                NULL, &dims[0], NULL);
        status = H5Dwrite(dataset, h5_atomic_type, dataspace, filespace,
                H5P_DEFAULT, data.origin());
        ++offset[data_rank];
    }

template<>
    void
    Hdf5_writer<MArray2d_ref >::append(MArray2d_ref & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(data_rank);
            for (int i = 0; i < data_rank; ++i) {
                data_dims.at(i) = data.shape()[i];
            }
            setup(data_dims, h5_atomic_typename<double > ());
        }
        ++size[data_rank];
        status = H5Dextend(dataset, &size[0]);

        filespace = H5Dget_space(dataset);
        have_filespace = true;
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0],
                NULL, &dims[0], NULL);
        status = H5Dwrite(dataset, h5_atomic_type, dataspace, filespace,
                H5P_DEFAULT, data.origin());
        ++offset[data_rank];
    }

template<>
    void
    Hdf5_writer<MArray3d_ref >::append(MArray3d_ref & data)
    {
        if (!have_setup) {
            std::vector<int > data_dims(data_rank);
            for (int i = 0; i < data_rank; ++i) {
                data_dims.at(i) = data.shape()[i];
            }
            setup(data_dims, h5_atomic_typename<double > ());
        }
        ++size[data_rank];
        status = H5Dextend(dataset, &size[0]);

        filespace = H5Dget_space(dataset);
        have_filespace = true;
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset[0],
                NULL, &dims[0], NULL);
        status = H5Dwrite(dataset, h5_atomic_type, dataspace, filespace,
                H5P_DEFAULT, data.origin());
        ++offset[data_rank];
    }
