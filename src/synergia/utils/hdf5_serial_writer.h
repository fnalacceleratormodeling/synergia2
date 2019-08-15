#ifndef HDF5_SERIAL_WRITER_H_
#define HDF5_SERIAL_WRITER_H_

#include <vector>
#include <string>

#include "hdf5.h"

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"
#include "synergia/utils/cereal.h"

class Hdf5_serial_writer
{
private:

    int data_rank;
    std::string name;
    std::vector<hsize_t> dims, max_dims, size, offset;

    hid_t file_ptr;
    Hdf5_handler dataset;
    Hdf5_handler atomic_type;

    bool have_setup;
    bool resume;
    size_t data_size;

private:

    void do_setup(std::vector<int> const& data_dims);
    void do_append(void* ptr);

    template<typename T>
    std::vector<int> pre_setup(T const& data)
    {
        data_rank = 0;
        atomic_type = hdf5_atomic_data_type<T>();
        data_size = sizeof(T);

        // dummy variable -- length really should
        // be 0, but that would not compile
        return std::vector<int>(1);
    }

    template<typename T>
    void* get_data_ptr(T const& data)
    { return (void*)&data; }

public:

    Hdf5_serial_writer(
            hid_t file_ptr,
            std::string const& name)
        : data_rank(0), name(name)
        , dims(), max_dims(), size(), offset()
        , file_ptr(file_ptr)
        , dataset()
        , atomic_type()
        , have_setup(false)
        , resume(false)
        , data_size(0)
    { }

    template<typename T>
    void append(T & data)
    {
        if (!have_setup) do_setup(pre_setup(data));
        do_append(get_data_ptr(data));
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar & CEREAL_NVP(name);
        ar & CEREAL_NVP(offset);
        //ar & CEREAL_NVP(file_ptr);
        ar & CEREAL_NVP(data_rank);
        ar & CEREAL_NVP(data_size);
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar & CEREAL_NVP(name);
        ar & CEREAL_NVP(offset);
        //ar & CEREAL_NVP(file_ptr);
        ar & CEREAL_NVP(data_rank);
        ar & CEREAL_NVP(data_size);
        resume = true;
        have_setup = false;
    }

};


template<>
inline std::vector<int> Hdf5_serial_writer::pre_setup<karray1d>(karray1d const& data)
{
    data_rank = 1;
    atomic_type = hdf5_atomic_data_type<double>();
    data_size = sizeof(double) * data.size();

    std::vector<int> data_dims(data_rank);
    for (int i=0; i<data_rank; ++i) data_dims[i] = data.extent(i);
    return data_dims;
}

template<>
inline void* Hdf5_serial_writer::get_data_ptr<karray1d>(karray1d const& data)
{
    return data.data();
}


#endif /* HDF5_SERIAL_WRITER_H_ */
