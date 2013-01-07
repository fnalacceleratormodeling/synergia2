#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"
#include <stdexcept>
#include <iostream>

#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif

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
            dims[i] = data_dims.at(i);
            max_dims[i] = data_dims.at(i);
            size[i] = data_dims.at(i);
            chunk_dims[i] = data_dims.at(i);
            offset[i] = 0;
        }
        max_dims[data_rank] = H5S_UNLIMITED;
        if (resume) {
            dataset = file_sptr->get_h5file().openDataSet(name.c_str());
            DataSpace dataspace = dataset.getSpace();
            int file_rank = dataspace.getSimpleExtentNdims();
            if (file_rank != data_rank + 1) {
                throw std::runtime_error(
                        "Hdf5_serial_writer::resumed data has wrong rank");
            }
            dataspace.getSimpleExtentDims(&dims[0], NULL);
            size[data_rank] = dims[data_rank];
            offset[data_rank] = dims[data_rank];
            dims[data_rank] = 1;
        } else {
            size[data_rank] = 0;
            offset[data_rank] = 0;
            dims[data_rank] = 1;
            if (data_size == 0) {
                throw std::runtime_error("Hdf5_serial_writer: zero data size encountered");
            }
            const size_t good_chunk_size = 8192; // pulled out of air
            if (data_size < good_chunk_size) {
                chunk_dims[data_rank] = good_chunk_size/data_size;
            } else {
                chunk_dims[data_rank] = 1;
            }
            DSetCreatPropList cparms;
            cparms.setChunk(data_rank + 1, &chunk_dims[0]);
            DataSpace dataspace(data_rank + 1, &dims[0], &max_dims[0]);
            dataset = file_sptr->get_h5file().createDataSet(name.c_str(),
                    atomic_type, dataspace, cparms);
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
        DataSpace dataspace(data_rank + 1, &dims[0], &max_dims[0]);
        ++size[data_rank];
        dataset.extend(&size[0]);

        DataSpace filespace = dataset.getSpace();
        filespace.selectHyperslab(H5S_SELECT_SET, &dims[0], &offset[0]);
        dataset.write(&data, atomic_type, dataspace, filespace);
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
    }
