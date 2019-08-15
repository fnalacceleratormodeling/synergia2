#ifndef HDF5_FILE_H_
#define HDF5_FILE_H_

#include "hdf5.h"

#include <string>
#include <memory>

#include "synergia/utils/hdf5_misc.h"
#include "synergia/utils/hdf5_writer.h"
#include "synergia/utils/hdf5_serial_writer.h"

#include "synergia/utils/cereal.h"
#include "synergia/utils/cereal_files.h"

#include "synergia/utils/commxx.h"


class Hdf5_file
{
public:

    enum Flag { truncate, read_write, read_only };
    enum Atomic_type { double_type, int_type };

private:

    std::string file_name;
    Hdf5_handler h5file;
    bool is_open;
    Flag current_flag;

    std::map<std::string, Hdf5_serial_writer> swriters;

    static unsigned int flag_to_h5_flags(Flag flag)
    {
        if (flag == Hdf5_file::truncate)   return H5F_ACC_TRUNC;
        if (flag == Hdf5_file::read_write) return H5F_ACC_RDWR;
        if (flag == Hdf5_file::read_only)  return H5F_ACC_RDONLY;

        return 0;
    }

public:

    Hdf5_file(std::string const& file_name, Flag flag);
    ~Hdf5_file() { close(); }

    void open(Flag flag);
    void close();
    void flush() const;

    std::vector<std::string> get_member_names();
    Atomic_type get_atomic_type(std::string const& name);
    std::vector<int > get_dims(std::string const& name);

    hid_t get_h5file()
    { return h5file.hid; }

    template<typename T>
    void write(T const& data, std::string const& name)
    { Hdf5_writer<T>(h5file.hid, name).write(data); }

    template<typename T>
    void write(T const* data, size_t len, std::string const & name)
    { Hdf5_writer<T>(h5file.hid, name).write(data, len); }

    template<typename T>
    void write_serial(T const& data, std::string const& name)
    { swriters.emplace(name, h5file.hid, name).first->second.append<T>(data); }

    template<typename T>
    T read(std::string const& name)
    {
        T retval;

        Hdf5_handler atomic_type = hdf5_atomic_data_type<T>();
        Hdf5_handler dataset = H5Dopen2(h5file.hid, name.c_str(), H5P_DEFAULT);

        std::vector<hsize_t > dims(1);
        dims.at(0) = 1;
        int data_rank = 0;

        Hdf5_handler dataspace = H5Screate_simple(data_rank, &dims[0], NULL);
        Hdf5_handler memspace = H5Screate_simple(data_rank, &dims[0], NULL);

        herr_t res = H5Dread(dataset.hid, atomic_type.hid, memspace.hid, 
                dataspace.hid, H5P_DEFAULT, &retval);

        if (res < 0) 
            throw Hdf5_exception("Error at reading Dataset " + name + " from HDF5 file");

        return retval;
    }


    template<class Archive>
    void save(Archive & ar, const unsigned int version) const;

    template<class Archive>
    void load(Archive & ar, const unsigned int version);

};

typedef std::shared_ptr<Hdf5_file > Hdf5_file_sptr; // syndoc:include


template<>
int* Hdf5_file::read<int*>(std::string const& name);

template<>
double* Hdf5_file::read<double*>(std::string const& name);

#if 0
template<>
MArray1d Hdf5_file::read<MArray1d >(std::string const& name);

template<>
MArray2d Hdf5_file::read<MArray2d >(std::string const& name);

template<>
MArray3d Hdf5_file::read<MArray3d >(std::string const& name);

template<>
MArray1i Hdf5_file::read<MArray1i >(std::string const& name);
#endif


#endif /* HDF5_FILE_H_ */
