#ifndef HDF5_FILE_H_
#define HDF5_FILE_H_
#include <string>
#include "hdf5.h"

#include "synergia/utils/hdf5_writer.h"

class Hdf5_file
{
private:
    hid_t h5file;
public:
    Hdf5_file(std::string const& file_name);
    template<typename T>
        void
        write(T const& data, std::string const& name)
        {
            Hdf5_writer<T > (h5file, name).write(data);
        }
    ;
    ~Hdf5_file();
};

#endif /* HDF5_FILE_H_ */
