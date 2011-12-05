#ifndef HDF5_FILE_H_
#define HDF5_FILE_H_
#include <string>
#include "H5Cpp.h"

#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"

class Hdf5_file
{
private:
    H5::H5File file;
public:
    Hdf5_file(std::string const& file_name, bool read_only = false);
    template<typename T>
        void
        write(T const& data, std::string const& name)
        {
            Hdf5_writer<T > (file, name).write(data);
        }
    ;
    template<typename T>
        T
        read(std::string const& name)
        {
            T retval;
            DataSet dataset = file.openDataSet(name.c_str());
            H5::DataType atomic_type = hdf5_atomic_data_type<T > ();
            std::vector<hsize_t > dims(1);
            dims.at(0) = 1;
            int data_rank = 0;
            DataSpace dataspace(data_rank, &dims[0]);
            DataSpace memspace(data_rank, &dims[0]);
            T * data_out = &retval;
            dataset.read(data_out, atomic_type, memspace, dataspace);
            return retval;
        }
    ;
    ~Hdf5_file();
};

#endif /* HDF5_FILE_H_ */
