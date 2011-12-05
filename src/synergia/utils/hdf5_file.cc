#include "hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"

Hdf5_file::Hdf5_file(std::string const& file_name, bool read_only) :
    file(file_name.c_str(), read_only ? H5F_ACC_RDONLY : H5F_ACC_TRUNC)
{
}

Hdf5_file::~Hdf5_file()
{
    file.close();
}

template<>
    MArray1d
    Hdf5_file::read<MArray1d >(std::string const& name)
    {
        DataSet dataset = file.openDataSet(name.c_str());
        H5::DataType atomic_type = hdf5_atomic_data_type<double > ();

        const int rank = 1;
        std::vector<hsize_t > dims(rank);
        DataSpace dataspace = dataset.getSpace();
        int file_rank = dataspace.getSimpleExtentNdims();
        if (file_rank != rank) {
            throw std::runtime_error(
                    "Hdf5_file::read<MArray1d>: data to read has wrong rank");
        }
        dataspace.getSimpleExtentDims(&dims[0], NULL);
        MArray1d retval(boost::extents[dims[0]]);

        DataSpace memspace(rank, &dims[0]);
        double * data_out = retval.origin();
        dataset.read(data_out, atomic_type, memspace, dataspace);
        return retval;
    }

template<>
    MArray2d
    Hdf5_file::read<MArray2d >(std::string const& name)
    {
        DataSet dataset = file.openDataSet(name.c_str());
        H5::DataType atomic_type = hdf5_atomic_data_type<double > ();

        const int rank = 2;
        std::vector<hsize_t > dims(rank);
        DataSpace dataspace = dataset.getSpace();
        int file_rank = dataspace.getSimpleExtentNdims();
        if (file_rank != rank) {
            throw std::runtime_error(
                    "Hdf5_file::read<MArray2d>: data to read has wrong rank");
        }
        dataspace.getSimpleExtentDims(&dims[0], NULL);
        MArray2d retval(boost::extents[dims[0]][dims[1]]);

        DataSpace memspace(rank, &dims[0]);
        double * data_out = retval.origin();
        dataset.read(data_out, atomic_type, memspace, dataspace);
        return retval;
    }

template<>
    MArray3d
    Hdf5_file::read<MArray3d >(std::string const& name)
    {
        DataSet dataset = file.openDataSet(name.c_str());
        H5::DataType atomic_type = hdf5_atomic_data_type<double > ();

        const int rank = 3;
        std::vector<hsize_t > dims(rank);
        DataSpace dataspace = dataset.getSpace();
        int file_rank = dataspace.getSimpleExtentNdims();
        if (file_rank != rank) {
            throw std::runtime_error(
                    "Hdf5_file::read<MArray3d>: data to read has wrong rank");
        }
        dataspace.getSimpleExtentDims(&dims[0], NULL);
        MArray3d retval(boost::extents[dims[0]][dims[1]][dims[2]]);

        DataSpace memspace(rank, &dims[0]);
        double * data_out = retval.origin();
        dataset.read(data_out, atomic_type, memspace, dataspace);
        return retval;
    }
