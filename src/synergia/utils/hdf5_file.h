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

    // set the type based on cmake config
#ifdef USE_PARALLEL_HDF5
    using use_parallel = std::true_type;
#else
    using use_parallel = std::false_type;
#endif

private:

    std::shared_ptr<Commxx> comm;
    std::string file_name;
    Hdf5_handler h5file;
    int root_rank;
    bool is_open;
    Flag current_flag;
    bool has_file;

    std::map<std::string, Hdf5_serial_writer> swriters;

    static unsigned int flag_to_h5_flags(Flag flag)
    {
        if (flag == Hdf5_file::truncate)   return H5F_ACC_TRUNC;
        if (flag == Hdf5_file::read_write) return H5F_ACC_RDWR;
        if (flag == Hdf5_file::read_only)  return H5F_ACC_RDONLY;

        return 0;
    }

public:

    Hdf5_file( std::string const& file_name, 
               Flag flag, 
               Commxx const& comm = Commxx() );

    ~Hdf5_file() { close(); }

    void open(Flag flag);
    void close();
    void flush() const;

#if 0
    std::vector<std::string> get_member_names();
    Atomic_type get_atomic_type(std::string const& name);
    std::vector<int > get_dims(std::string const& name);
#endif

    hid_t get_h5file()
    { return h5file.hid; }

    int master_rank() const
    { return root_rank; }

    // gather on the first dimension. all other dimensions must be of the same extents
    // calling from 4 ranks:
    //   write_collective("ds", pz) -> "ds" : [pz, pz, ...]
    //   write_collective("part", part[0:1][0:6]) -> "part" : part[0:3][0:6]
    template<typename T>
    void write_collective(std::string const& name, T const& data)
    { write(name, data, true); }

    // no gather, only the root rank will execute the write
    // calling from 4 ranks:
    //   write_single("ds", pz) -> "ds" : pz
    //   write_single("part", part[0:1][0:6]) -> "part" : part[0:1][0:6]
    template<typename T>
    void write_single(std::string const& name, T const& data)
    { write(name, data, false); }

    template<typename T>
    void write(std::string const& name, T const& data, bool collective = false)
    { Hdf5_writer::write(h5file, name, data, collective, *comm, root_rank); }

    template<typename T>
    void write(std::string const & name, T const* data, size_t len, bool collective = false)
    { Hdf5_writer::write(h5file, name, data, len, collective, *comm, root_rank); }


    // same as write_single(), except this will do append instead of overwrite
    template<typename T>
    void append_single(std::string const& name, T const& data)
    { }

    // same as write(), except this will do append instead of overwrite
    template<typename T>
    void append(std::string const& name, T const& data)
    { 
        auto w = swriters.find(name);

        if (w == swriters.end()) 
            w = swriters.emplace(name, Hdf5_serial_writer(h5file.hid, name)).first;

        w->second.append(data);
    }

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

    template<typename T>
    void read_array(std::string const& name, T * const ptr);

private:

    friend class cereal::access;

    Hdf5_file() 
    : file_name(), h5file(), is_open(false), current_flag(Hdf5_file::read_only) 
    { }

    template<class Archive>
    void save(Archive & ar) const
    {
        ar(CEREAL_NVP(comm));
        ar(CEREAL_NVP(file_name));
        ar(CEREAL_NVP(root_rank));
        ar(CEREAL_NVP(is_open));
        ar(CEREAL_NVP(current_flag));

        if (is_open)
        {
            flush();
            // TODO: copy_to_serialization_directory(file_name);
        }
    }

    template<class Archive>
    void load(Archive & ar)
    {
        ar(CEREAL_NVP(comm));
        ar(CEREAL_NVP(file_name));
        ar(CEREAL_NVP(root_rank));
        ar(CEREAL_NVP(is_open));
        ar(CEREAL_NVP(current_flag));

        if (is_open)
        {
            // TODO: copy_from_serialization_directory(file_name)
            is_open = false;
            open(current_flag);
        }
    }
};


template<>
int* Hdf5_file::read<int*>(std::string const& name);

template<>
double* Hdf5_file::read<double*>(std::string const& name);

template<>
karray2d Hdf5_file::read<karray2d>(std::string const& name);

template<>
karray2d_row Hdf5_file::read<karray2d_row>(std::string const& name);

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
