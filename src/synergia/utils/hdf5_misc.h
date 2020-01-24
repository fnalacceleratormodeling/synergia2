#ifndef HDF5_MISC_H_
#define HDF5_MISC_H_

#include "hdf5.h"

#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <exception>

namespace storage_order
{
    constexpr int c = 0;
    constexpr int fortran = 1;

    constexpr int right = c;
    constexpr int left = fortran;

    constexpr int row = c;
    constexpr int col = fortran;

    constexpr int hdf5_default = c;
}

// hdf5 exceptions
struct Hdf5_exception : public std::exception
{
    Hdf5_exception(std::string const & msg = "")
        : hdf5_msg(), user_msg(msg)
    { 
        std::stringstream buf;

        //H5Eprint(H5E_DEFAULT, stderr);
        H5Ewalk(H5E_DEFAULT, H5E_WALK_UPWARD, &Hdf5_exception::err_walk, &buf);

        hdf5_msg = buf.str();
    }

    ~Hdf5_exception() throw() { }

    virtual const char * what() const throw()
    {
        std::string res = std::string("\n")
            + "===================================\n" 
            + "USER MESSAGE:\n  " + user_msg + "\n\n"
            + "HDF5 MESSAGE:\n"   + hdf5_msg
            + "===================================\n";
        return res.c_str();
    }

private:

    static herr_t err_walk(unsigned int n, H5E_error_t const *err_desc, void *client_data)
    {
        std::stringstream *ss = (std::stringstream *)client_data;

        /* Check arguments */
        if (!client_data) return 0;

        /* Get descriptions for the major and minor error numbers */
        const char *maj_str = H5Eget_major(err_desc->maj_num);
        const char *min_str = H5Eget_minor(err_desc->min_num);

        /* Print error message */
        (*ss) << "  #" << std::setfill('0') << std::setw(4) << n
              << ": " << err_desc->file_name << " line " << err_desc->line 
              << " in " << err_desc->func_name << "(): " << err_desc->desc << "\n"
              << "    major(" << std::setfill('0') << std::setw(3) << (int)err_desc->maj_num 
              << "): " << maj_str << "\n"
              << "    minor(" << std::setfill('0') << std::setw(3) << (int)err_desc->min_num 
              << "): " << min_str << "\n"
            ;

        return 0;
    }

    std::string hdf5_msg;
    std::string user_msg;
};

// handles the closing of resources in the RAII way
struct Hdf5_handler
{
    Hdf5_handler(hid_t handler = 0) : hid(handler)
    { if (hid < 0) throw Hdf5_exception("Bad HDF5 Handler"); }

    Hdf5_handler(Hdf5_handler&& o) noexcept : hid(o.hid)
    { o.hid = 0; }

    Hdf5_handler & operator= (hid_t handler)
    {
        if (handler < 0) throw Hdf5_exception("Bad HDF5 Handler");
        hid = handler; return *this;
    }

    // no copy and copy assignment
    Hdf5_handler(Hdf5_handler const&) = delete;
    Hdf5_handler& operator= (Hdf5_handler const&) = delete;

    ~Hdf5_handler()
    { close(); }

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
        default: break;
        }

        hid = 0;
    }

    operator hid_t() const
    { return hid; }

    hid_t hid;
};


// The generic (T) version of h5_atomic_data_type is undefined.
// Only versions with specializations will compile.
template<typename T>
inline hid_t hdf5_atomic_data_type();

template<>
inline hid_t hdf5_atomic_data_type<uint8_t>()
{ return H5Tcopy(H5T_NATIVE_UCHAR); }

template<>
inline hid_t hdf5_atomic_data_type<int>()
{ return H5Tcopy(H5T_NATIVE_INT); }

template<>
inline hid_t hdf5_atomic_data_type<double>()
{ return H5Tcopy(H5T_NATIVE_DOUBLE); }



#endif /* HDF5_MISC_H_ */
