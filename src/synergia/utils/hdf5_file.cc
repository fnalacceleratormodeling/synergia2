#include "hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"

// flag_to_h5_flags is a local function
unsigned int
flag_to_h5_flags(Hdf5_file::Flag flag)
{
    unsigned int retval;
    if (flag == Hdf5_file::truncate) {
        retval = H5F_ACC_TRUNC;
    } else if (flag == Hdf5_file::read_write) {
        retval = H5F_ACC_RDWR;
    } else if (flag == Hdf5_file::read_only) {
        retval = H5F_ACC_RDONLY;
    } else {
        retval = 0;
    }
    return retval;
}

Hdf5_file::Hdf5_file(std::string const& file_name, Flag flag) :
    file_name(file_name), h5file_ptr(0), is_open(false)
{
    open(flag);
}

void
Hdf5_file::open(Flag flag)
{
    if (is_open) {
        close();
    }
    unsigned int h5_flag;
    if (flag == Hdf5_file::truncate) {
        h5_flag = H5F_ACC_TRUNC;
    } else if (flag == Hdf5_file::read_write) {
        h5_flag = H5F_ACC_RDWR;
    } else if (flag == Hdf5_file::read_only) {
        h5_flag = H5F_ACC_RDONLY;
    }
    h5file_ptr = new H5::H5File(file_name.c_str(), h5_flag);
    is_open = true;
}

void
Hdf5_file::close()
{
    if (is_open) {
        h5file_ptr->close();
        delete h5file_ptr;
        is_open = false;
    }
}

Hdf5_file::~Hdf5_file()
{
    close();
}

template<>
    MArray1d
    Hdf5_file::read<MArray1d >(std::string const& name)
    {
        DataSet dataset = h5file_ptr->openDataSet(name.c_str());
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
        DataSet dataset = h5file_ptr->openDataSet(name.c_str());
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
        DataSet dataset = h5file_ptr->openDataSet(name.c_str());
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
