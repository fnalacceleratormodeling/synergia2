#ifndef HDF5_FILE_H_
#define HDF5_FILE_H_
#include <string>
#include "H5Cpp.h"

#include "synergia/utils/hdf5_writer.h"

class Hdf5_file
{
public:
    enum Flag
    {
        truncate, read_write, read_only
    };
private:
    std::string file_name;
    H5::H5File * h5file_ptr;
    bool is_open;
public:
    Hdf5_file(std::string const& file_name, Flag flag);
    void
    open(Flag flag);
    void
    close();
    H5::H5File &
    get_h5file();
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
