#include "hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"

Hdf5_file::Hdf5_file( 
        std::string const& file_name, 
        Flag flag, 
        Commxx const& c )
    : comm(std::make_shared<Commxx>(c))
    , file_name(file_name)
    , h5file()
    , root_rank(c.size()-1)
    , is_open(false)
    , current_flag(flag)
#ifdef USE_PARALLEL_HDF5
    , has_file(true)
#else
    , has_file(c.rank() == root_rank)
#endif
{
    // turn off the automatic error printing
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    // open file
    open(flag);
}

Hdf5_file::Hdf5_file( 
        std::string const& file_name, 
        Flag flag, 
        std::shared_ptr<Commxx> const& c )
    : comm(c)
    , file_name(file_name)
    , h5file()
    , root_rank(c->size()-1)
    , is_open(false)
    , current_flag(flag)
#ifdef USE_PARALLEL_HDF5
    , has_file(true)
#else
    , has_file(c->rank() == root_rank)
#endif
{
    // turn off the automatic error printing
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    // open file
    open(flag);
}

void Hdf5_file::open(Flag flag)
{
    // open file on this rank?
    if (!has_file) return;

    // already opened
    if (is_open) close();

    int attempts = 0;
    bool fail = true;

    while ((attempts < 5) && fail)
    {
        try
        {
            Hdf5_handler plist_id = H5Pcreate(H5P_FILE_ACCESS);

#ifdef USE_PARALLEL_HDF5
            MPI_Info info = MPI_INFO_NULL;
            H5Pset_fapl_mpio(plist_id, *comm, info);
#endif

            if (flag == Hdf5_file::truncate)
            {
                // create
                h5file = H5Fcreate(file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
                current_flag = Hdf5_file::read_write;
            }
            else
            {
                // open
                h5file = H5Fopen(file_name.c_str(), flag_to_h5_flags(flag), plist_id);
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

#if 0
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
#endif


#if 0
template<>
void Hdf5_file::read_array<double>(std::string const& name, double * const ptr)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<double>();

    const int rank = 1;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read_array<double>: data to read has wrong rank");
    }

    // dims
    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception("error getting dims from dataspace");

    // memspace
    Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);

    // read
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, ptr);
    if (res < 0) throw Hdf5_exception("error reading from hdf5 file");
}

template<>
void Hdf5_file::read_array<uint8_t>(std::string const& name, uint8_t* const ptr)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<uint8_t>();

    const int rank = 1;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read_array<uint8_t>: data to read has wrong rank");
    }

    // dims
    std::vector<hsize_t> dims(rank);
    herr_t res = H5Sget_simple_extent_dims(dataspace, &dims[0], NULL);
    if (res < 0) throw Hdf5_exception("error getting dims from dataspace");

    // memspace
    Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);

    // read
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, ptr);
    if (res < 0) throw Hdf5_exception("error reading from hdf5 file");
}

template<>
karray2d Hdf5_file::read<karray2d>(std::string const& name)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<double>();

    const int rank = 2;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read<karray2d>: data to read has wrong rank");
    }

    hsize_t dims[rank];
    herr_t res = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    if (res < 0) throw Hdf5_exception();

    int storage_order = storage_order::hdf5_default;
    try { storage_order = read<int>(name + "_storage_order"); } catch (Hdf5_exception e) { }

    if (storage_order != storage_order::col)
    {
        throw std::runtime_error(
                "Hdf5_file::read<karray2d>: data to read has wrong storage order");
    }

    karray2d retval(name, dims[0], dims[1]);

    Hdf5_handler memspace = H5Screate_simple(rank, dims, NULL);
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval.data());

    if (res < 0) throw Hdf5_exception();

    return retval;
}

template<>
karray2d_row Hdf5_file::read<karray2d_row>(std::string const& name)
{
    Hdf5_handler dataset     = H5Dopen(h5file, name.c_str(), H5P_DEFAULT);
    Hdf5_handler dataspace   = H5Dget_space(dataset);
    Hdf5_handler atomic_type = hdf5_atomic_data_type<double>();

    const int rank = 2;
    int file_rank = H5Sget_simple_extent_ndims(dataspace);

    if (file_rank != rank)
    {
        throw std::runtime_error(
                "Hdf5_file::read<karray2d_row>: data to read has wrong rank");
    }

    hsize_t dims[rank];
    herr_t res = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    if (res < 0) throw Hdf5_exception();

    int storage_order = storage_order::hdf5_default;
    try { storage_order = read<int>(name + "_storage_order"); } catch (Hdf5_exception e) { }

    if (storage_order != storage_order::row)
    {
        throw std::runtime_error(
                "Hdf5_file::read<karray2d_row>: data to read has wrong storage order");
    }

    karray2d_row retval(name, dims[0], dims[1]);

    Hdf5_handler memspace = H5Screate_simple(rank, dims, NULL);
    res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval.data());

    if (res < 0) throw Hdf5_exception();

    return retval;
}
#endif



