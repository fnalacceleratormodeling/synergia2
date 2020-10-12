#ifndef HDF5_MISC_H_
#define HDF5_MISC_H_

#include "hdf5.h"

#include <stdexcept>
#include <string>

// hdf5 exceptions
struct Hdf5_exception : public std::runtime_error
{
  explicit Hdf5_exception(std::string const& msg) : std::runtime_error(msg)
    {
      H5Eclear(H5E_DEFAULT);
    }

    explicit Hdf5_exception(char const* msg = "")
        : std::runtime_error(msg)
    {
        H5Eclear(H5E_DEFAULT);
    }
};

// handles the closing of resources in the RAII way
struct Hdf5_handler
{
    Hdf5_handler(hid_t handler = 0) : hid(handler)
    { 
        if (hid < 0) throw Hdf5_exception("Bad HDF5 Handler");
    }

    ~Hdf5_handler()
    {
        close();
    }

    void close()
    {
        switch(H5Iget_type(hid))
        {
        case H5I_FILE:      H5Fclose(hid); break;
        case H5I_GROUP:     H5Gclose(hid); break;
        case H5I_DATATYPE:  H5Tclose(hid); break;
        case H5I_DATASPACE: H5Sclose(hid); break;
        case H5I_DATASET:   H5Dclose(hid); break;
        case H5I_ATTR:      H5Aclose(hid); break;
        case H5I_GENPROP_LST: H5Pclose(hid); break;
        default:            break;
        }

        hid = 0;
    }

    Hdf5_handler & operator= (hid_t handler)
    {
        hid = handler; 
        if (hid < 0) throw Hdf5_exception("Bad HDF5 Handler");

        return *this;
    }

    operator hid_t () const
    {
        return hid;
    }

    hid_t hid;

private:

    // disable copy and assignment
    // Hdf5_handler(Hdf5_handler const &) { }
  // TODO: use C++11 delete'd function
    Hdf5_handler & operator= (Hdf5_handler const &);
};


// The generic (T) version of h5_atomic_data_type is undefined.
// Only versions with specializations will compile.
template<typename T>
    inline hid_t
    hdf5_atomic_data_type();

template<>
    inline hid_t
    hdf5_atomic_data_type<int > ()
    {
        return H5Tcopy(H5T_NATIVE_INT);
    }

template<>
    inline hid_t
    hdf5_atomic_data_type<double > ()
    {
        return H5Tcopy(H5T_NATIVE_DOUBLE);
    }



#endif /* HDF5_MISC_H_ */
