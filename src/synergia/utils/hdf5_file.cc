#include "hdf5_file.h"

Hdf5_file::Hdf5_file(std::string const& file_name, bool read_only) :
    file(file_name.c_str(), read_only ? H5F_ACC_RDONLY : H5F_ACC_TRUNC)
{
}

Hdf5_file::~Hdf5_file()
{
    file.close();
}
