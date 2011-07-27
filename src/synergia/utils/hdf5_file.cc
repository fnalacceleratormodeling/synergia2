#include "hdf5_file.h"

Hdf5_file::Hdf5_file(std::string const& file_name) :
    file(file_name.c_str(), H5F_ACC_TRUNC)
{
}

Hdf5_file::~Hdf5_file()
{
    file.close();
}
