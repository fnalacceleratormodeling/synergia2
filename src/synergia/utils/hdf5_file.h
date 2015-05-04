#ifndef HDF5_FILE_H_
#define HDF5_FILE_H_
#include <string>
#include "H5Cpp.h"
#include <boost/shared_ptr.hpp>

#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/serialization.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/commxx.h"

class Hdf5_file
{
public:
    enum Flag
    {
        truncate, read_write, read_only
    };
    enum Atomic_type
    {
        double_type, int_type
    };

private:
    std::string file_name;
    H5::H5File * h5file_ptr;
    bool is_open;
    Flag current_flag;
    unsigned int
    flag_to_h5_flags(Flag flag)
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
public:
    Hdf5_file(std::string const& file_name, Flag flag);
    // Default constructor for serialization use only
    Hdf5_file();
    void
    open(Flag flag);
    void
    close();
    void
    flush() const;
    std::vector<std::string>
    get_member_names();
    Atomic_type
    get_atomic_type(std::string const& name);
    std::vector<hsize_t>
    get_dims(std::string const& name);
    H5::H5File &
    get_h5file();
    template<typename T>
        void
        write(T const& data, std::string const& name);
    template<typename T>
        T
        read(std::string const& name);
    template<class Archive>
        void
        save(Archive & ar, const unsigned int version) const;
    template<class Archive>
        void
        load(Archive & ar, const unsigned int version);
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    ~Hdf5_file();
};

typedef boost::shared_ptr<Hdf5_file > Hdf5_file_sptr; // syndoc:include

#include "synergia/utils/hdf5_file.tcc"

#endif /* HDF5_FILE_H_ */
