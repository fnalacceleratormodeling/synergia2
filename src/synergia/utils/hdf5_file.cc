#include "hdf5_file.h"

Hdf5_file::Hdf5_file(std::string const& file_name)
{
    h5file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
            H5P_DEFAULT);
}

Hdf5_file::~Hdf5_file()
{
    H5Fclose(h5file);
}
