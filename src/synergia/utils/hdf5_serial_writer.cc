#include "hdf5_serial_writer.h"

template<>
    Hdf5_serial_writer<MArray1d_ref >::Hdf5_serial_writer(H5File & file,
            std::string const& name) :
        data_rank(MArray1d_ref::dimensionality), name(name), file(file),
                have_setup(false)
    {
    }

template<>
    Hdf5_serial_writer<MArray2d_ref >::Hdf5_serial_writer(H5File & file,
            std::string const& name) :
        data_rank(MArray2d_ref::dimensionality), name(name), file(file),
                have_setup(false)
    {
    }

template<>
    Hdf5_serial_writer<MArray3d_ref >::Hdf5_serial_writer(H5File & file,
            std::string const& name) :
        data_rank(MArray3d_ref::dimensionality), name(name), file(file),
                have_setup(false)
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
            setup(data_dims, hdf5_atomic_data_type<double > ());
        }
        DataSpace dataspace(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];
        dataset.extend(&size[0]);

        DataSpace filespace = dataset.getSpace();
        have_filespace = true;
        filespace.selectHyperslab(H5S_SELECT_SET, &dims[0], &offset[0]);
        dataset.write(data.origin(), h5_atomic_type, dataspace, filespace);
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
            setup(data_dims, hdf5_atomic_data_type<double > ());
        }
        DataSpace dataspace(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];
        dataset.extend(&size[0]);

        DataSpace filespace = dataset.getSpace();
        have_filespace = true;
        filespace.selectHyperslab(H5S_SELECT_SET, &dims[0], &offset[0]);
        dataset.write(data.origin(), h5_atomic_type, dataspace, filespace);
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
            setup(data_dims, hdf5_atomic_data_type<double > ());
        }
        DataSpace dataspace(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];
        dataset.extend(&size[0]);

        DataSpace filespace = dataset.getSpace();
        have_filespace = true;
        filespace.selectHyperslab(H5S_SELECT_SET, &dims[0], &offset[0]);
        dataset.write(data.origin(), h5_atomic_type, dataspace, filespace);
        ++offset[data_rank];
    }
