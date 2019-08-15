#include "hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"
//#include <boost/align/aligned_alloc.hpp>

Hdf5_file::Hdf5_file(std::string const& file_name, Flag flag) 
    : file_name(file_name)
    , h5file()
    , is_open(false)
{
    // turn off error printing
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    open(flag);
    current_flag = flag;
}

void Hdf5_file::open(Flag flag)
{
    if (is_open) close();

    int attempts = 0;
    bool fail = true;

    while ((attempts < 5) && fail)
    {
        try
        {
            if (flag == Hdf5_file::truncate)
            {
                // create
                h5file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            }
            else
            {
                // open
                h5file = H5Fopen(file_name.c_str(), flag_to_h5_flags(flag), H5P_DEFAULT);
            }

            fail = false;
        }
        catch(Hdf5_exception & e)
        {
            ++attempts;
            fail = true;

            std::cout << e.what() << "\n";
            std::cout << "caught hdf5 open file error, attempts number="
                << attempts << " on rank=" << Commxx().rank() << std::endl;

            sleep(3);
        }
    }

    is_open = true;
}

void Hdf5_file::close()
{
    if (is_open) {
        h5file.close();
        is_open = false;
    }
}

void
Hdf5_file::flush() const
{
    if (is_open) H5Fflush(h5file.hid, H5F_SCOPE_GLOBAL);
}

extern "C" herr_t get_member_names_callback(hid_t group, const char *name,
                                            const H5L_info_t *info,
                                            void *op_data)
{
    ((std::vector<std::string> *)op_data)->push_back(name);
    return 0;
}

std::vector<std::string>
Hdf5_file::get_member_names()
{
    std::vector<std::string> member_names;

    Hdf5_handler root_group = H5Gopen(h5file, "/", H5P_DEFAULT);
    herr_t status = H5Literate(root_group, H5_INDEX_NAME,
                               H5_ITER_NATIVE, NULL,
                               &get_member_names_callback,
                               (void *) &member_names);

    if (status < 0) throw Hdf5_exception("error at get hdf5 member names");

    return member_names;
}

Hdf5_file::Atomic_type
Hdf5_file::get_atomic_type(std::string const& name)
{
    Hdf5_handler dataset  = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler datatype = H5Dget_type(dataset);
    H5T_class_t the_class = H5Tget_class(datatype);

    Hdf5_file::Atomic_type retval;

    if (the_class == H5T_FLOAT)
    {
        retval = Hdf5_file::double_type;
    }
    else if (the_class == H5T_INTEGER)
    {
        retval = Hdf5_file::int_type;
    }
    else
    {
        throw std::runtime_error("Hdf5_file::get_atomic_type: type not handled");
    }

    return retval;
}

std::vector<int>
Hdf5_file::get_dims(std::string const& name)
{
    Hdf5_handler dataset   = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace = H5Dget_space(dataset);

    int rank = H5Sget_simple_extent_ndims(dataspace);
    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
    if (res<0) throw Hdf5_exception("error at getting dims");

    std::vector<int > retval(rank);
    std::copy(dims.begin(), dims.end(), retval.begin());

    return retval;
}


template<>
int* Hdf5_file::read<int*>(std::string const& name)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<int>();

    const int rank = 1;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read<int *>: data to read has wrong rank");
    }

    // dims
    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception("error getting dims from dataspace");

    // memspace
    Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);

    // allocate
    int * retval = new int[dims[0]];

    // read
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval);
    if (res < 0) throw Hdf5_exception("error reading from hdf5 file");

    return retval;
}


template<>
double* Hdf5_file::read<double*>(std::string const& name)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<double>();

    const int rank = 1;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read<double *>: data to read has wrong rank");
    }

    // dims
    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception("error getting dims from dataspace");

    // memspace
    Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);

    // allocate
    double * retval = new double[dims[0] * sizeof(double)];

    // read
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval);
    if (res < 0) throw Hdf5_exception("error reading from hdf5 file");

    return retval;
}

#if 0
template<>
MArray1d Hdf5_file::read<MArray1d >(std::string const& name)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<double>();

    const int rank = 1;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read<MArray1d>: data to read has wrong rank");
    }

    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);

    if (res < 0) throw Hdf5_exception("error getting dims from dataspace");

    MArray1d retval(boost::extents[dims[0]]);

    Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval.origin());

    if (res < 0) throw Hdf5_exception("error reading from hdf5 file");

    return retval;
}

template<>
MArray2d
Hdf5_file::read<MArray2d >(std::string const& name)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<double>();

    const int rank = 2;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read<MArray2d>: data to read has wrong rank");
    }

    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception();

    int storage_order = Hdf5_writer<MArray2d>::c_storage_order;
    try { storage_order = read<int>(name + "_storage_order"); } catch (Hdf5_exception e) { }

    MArray2d retval = (storage_order == Hdf5_writer<MArray2d >::c_storage_order) ?  
        MArray2d(boost::extents[dims[0]][dims[1]], boost::c_storage_order()) :
        MArray2d(boost::extents[dims[0]][dims[1]], boost::fortran_storage_order());

    Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval.origin());

    if (res < 0) throw Hdf5_exception();

    return retval;
}

template<>
MArray3d
Hdf5_file::read<MArray3d >(std::string const& name)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<double>();

    const int rank = 3;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read<MArray3d>: data to read has wrong rank");
    }

    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception();

    int storage_order = Hdf5_writer<MArray3d>::c_storage_order;
    try { storage_order = read<int>(name + "_storage_order"); } catch (Hdf5_exception e) { } 

    MArray3d retval = (storage_order == Hdf5_writer<MArray3d >::c_storage_order) ?  
        MArray3d(boost::extents[dims[0]][dims[1]][dims[2]], boost::c_storage_order()) :
        MArray3d(boost::extents[dims[0]][dims[1]][dims[2]], boost::fortran_storage_order());

    Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval.origin());

    if (res < 0) throw Hdf5_exception();

    return retval;
}

template<>
MArray1i
Hdf5_file::read<MArray1i >(std::string const& name)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<int>();

    const int rank = 1;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read<MArray1i>: data to read has wrong rank");
    }

    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception();

    MArray1i retval(boost::extents[dims[0]]);

    Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval.origin());

    if (res < 0) throw Hdf5_exception();

    return retval;
}
#endif

#if 0
template<class Archive>
void
Hdf5_file::save(Archive & ar, const unsigned int version) const
{
    ar << CEREAL_NVP(file_name)
       << CEREAL_NVP(is_open)
       << CEREAL_NVP(current_flag);
    if (is_open) {
        flush();
        copy_to_serialization_directory(file_name);
    }
}

template<class Archive>
void
Hdf5_file::load(Archive & ar, const unsigned int version)
{
    ar >> CEREAL_NVP(file_name)
            >> CEREAL_NVP(is_open)
            >> CEREAL_NVP(current_flag);
    if (is_open) {
        copy_from_serialization_directory(file_name);
        Flag flag;
        if (current_flag == read_only) {
            flag = read_only;
        } else {
            flag = read_write;
        }
        is_open = false; // will be changed to true by open
        open(flag);
    }
}

template
void
Hdf5_file::save<cereal::BinaryOutputArchive>(
        cereal::BinaryOutputArchive & ar, const unsigned int version) const;

template
void
Hdf5_file::save<cereal::XMLOutputArchive>(
        cereal::XMLOutputArchive & ar, const unsigned int version) const;

template
void
Hdf5_file::load<cereal::BinaryInputArchive>(
        cereal::BinaryInputArchive & ar, const unsigned int version);

template
void
Hdf5_file::load<cereal::XMLInputArchive>(
        cereal::XMLInputArchive & ar, const unsigned int version);
#endif

