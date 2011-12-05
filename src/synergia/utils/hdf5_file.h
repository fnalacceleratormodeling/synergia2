#ifndef HDF5_FILE_H_
#define HDF5_FILE_H_
#include <string>
#include "H5Cpp.h"

#include "synergia/utils/hdf5_writer.h"

class Hdf5_file
{
private:
    H5::H5File file;
public:
    Hdf5_file(std::string const& file_name, bool read_only = false);
    template<typename T>
        void
        write(T const& data, std::string const& name);
    template<typename T>
        T
        read(std::string const& name);
    ~Hdf5_file();
};

#include "synergia/utils/hdf5_file.tcc"

#endif /* HDF5_FILE_H_ */
