#include "hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/hdf5_misc.h"

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

Hdf5_file::Hdf5_file() : h5file()
{
}

void
Hdf5_file::open(Flag flag)
{
    if (is_open) {
        close();
    }

#if 0
    int attempts=0;
    bool fail=true;
    while ((attempts<5) && fail){
        try{
            h5file_ptr = new H5::H5File(file_name.c_str(), flag_to_h5_flags(flag));
            fail=false;
        }
        catch (H5::FileIException& e) {
            ++attempts;
            fail=true;
            std::cout << "caught hdf5 open file error, attempts number="
                <<attempts<<" on rank="<<Commxx().get_rank()<<std::endl;
            if (h5file_ptr) delete h5file_ptr;
            sleep(3);
        }
    }
#endif

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
        }
        catch(Hdf5_exception & e)
        {
            ++attempts;
            fail = true;
            std::cout << "caught hdf5 open file error, attempts number="
                << attempts << " on rank=" << Commxx().get_rank() << std::endl;
            sleep(3);
        }
    }

    is_open = true;
}

void
Hdf5_file::close()
{
    if (is_open) {
#if 0
        h5file_ptr->close();
        delete h5file_ptr;
#endif
        h5file.close();
        is_open = false;
    }
}

void
Hdf5_file::flush() const
{
    //h5file_ptr->flush(H5F_SCOPE_GLOBAL);
    H5Fflush(h5file.hid, H5F_SCOPE_GLOBAL);
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
#if 0
    Group *root_group_ptr = new Group (h5file_ptr->openGroup("/"));
    herr_t status = H5Literate(root_group_ptr->getId(), H5_INDEX_NAME,
                               H5_ITER_NATIVE, NULL,
                               &get_member_names_callback,
                               (void *) &member_names);
    delete root_group_ptr;
#endif

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
#if 0
    DataSet dataset = h5file_ptr->openDataSet(name.c_str());
    DataType datatype = dataset.getDataType();
    H5T_class_t the_class = datatype.getClass();

    Hdf5_file::Atomic_type retval;
    if(the_class == H5T_FLOAT) {
        retval = Hdf5_file::double_type;
    } else if(the_class == H5T_INTEGER) {
        retval = Hdf5_file::int_type;
    } else {
        throw std::runtime_error("Hdf5_file::get_atomic_type: type not handled");
    }
    return retval;
#endif

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

std::vector<int >
Hdf5_file::get_dims(std::string const& name)
{
#if 0
    DataSet dataset = h5file_ptr->openDataSet(name.c_str());
    DataSpace dataspace = dataset.getSpace();
    int rank = dataspace.getSimpleExtentNdims();
    std::vector<hsize_t > dims(rank);
    dataspace.getSimpleExtentDims(&dims[0], NULL);

    std::vector<int > retval(rank);
    std::copy(dims.begin(), dims.end(), retval.begin());
    return retval;
#endif

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

hid_t
Hdf5_file::get_h5file()
{
    //return *h5file_ptr;
    return h5file.hid;
}

Hdf5_file::~Hdf5_file()
{
    close();
}

template<>
    MArray1d
    Hdf5_file::read<MArray1d >(std::string const& name)
    {
#if 0
        DataSet dataset = h5file_ptr->openDataSet(name.c_str());
        H5::DataType atomic_type = hdf5_atomic_data_type<double > ();

        const int rank = 1;
        std::vector<hsize_t > dims(rank);
        DataSpace dataspace = dataset.getSpace();
        int file_rank = dataspace.getSimpleExtentNdims();
        if (file_rank != rank) {
            throw std::runtime_error(
                    "Hdf5_file::read<MArray1d>: data to read has wrong rank");
        }
        dataspace.getSimpleExtentDims(&dims[0], NULL);
        MArray1d retval(boost::extents[dims[0]]);

        DataSpace memspace(rank, &dims[0]);
        double * data_out = retval.origin();
        dataset.read(data_out, atomic_type, memspace, dataspace);
        return retval;
#endif

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
#if 0
        DataSet dataset = h5file_ptr->openDataSet(name.c_str());
        H5::DataType atomic_type = hdf5_atomic_data_type<double > ();

        const int rank = 2;
        std::vector<hsize_t > dims(rank);
        DataSpace dataspace = dataset.getSpace();
        int file_rank = dataspace.getSimpleExtentNdims();
        if (file_rank != rank) {
            throw std::runtime_error(
                    "Hdf5_file::read<MArray2d>: data to read has wrong rank");
        }
        dataspace.getSimpleExtentDims(&dims[0], NULL);
        MArray2d retval(boost::extents[dims[0]][dims[1]]);

        DataSpace memspace(rank, &dims[0]);
        double * data_out = retval.origin();
        dataset.read(data_out, atomic_type, memspace, dataspace);
        return retval;
#endif

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

        MArray2d retval(boost::extents[dims[0]][dims[1]]);

        Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);
        res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval.origin());

        if (res < 0) throw Hdf5_exception();

        return retval;
    }

template<>
    MArray3d
    Hdf5_file::read<MArray3d >(std::string const& name)
    {
#if 0
        DataSet dataset = h5file_ptr->openDataSet(name.c_str());
        H5::DataType atomic_type = hdf5_atomic_data_type<double > ();

        const int rank = 3;
        std::vector<hsize_t > dims(rank);
        DataSpace dataspace = dataset.getSpace();
        int file_rank = dataspace.getSimpleExtentNdims();
        if (file_rank != rank) {
            throw std::runtime_error(
                    "Hdf5_file::read<MArray3d>: data to read has wrong rank");
        }
        dataspace.getSimpleExtentDims(&dims[0], NULL);
        MArray3d retval(boost::extents[dims[0]][dims[1]][dims[2]]);

        DataSpace memspace(rank, &dims[0]);
        double * data_out = retval.origin();
        dataset.read(data_out, atomic_type, memspace, dataspace);
        return retval;
#endif

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

        MArray3d retval(boost::extents[dims[0]][dims[1]][dims[2]]);

        Hdf5_handler memspace = H5Screate_simple(rank, &dims[0], NULL);
        res = H5Dread(dataset, atomic_type, memspace, dataspace, H5P_DEFAULT, retval.origin());

        if (res < 0) throw Hdf5_exception();

        return retval;
    }

template<>
    MArray1i
    Hdf5_file::read<MArray1i >(std::string const& name)
    {
#if 0
        DataSet dataset = h5file_ptr->openDataSet(name.c_str());
        H5::DataType atomic_type = hdf5_atomic_data_type<int > ();

        const int rank = 1;
        std::vector<hsize_t > dims(rank);
        DataSpace dataspace = dataset.getSpace();
        int file_rank = dataspace.getSimpleExtentNdims();
        if (file_rank != rank) {
            throw std::runtime_error(
                    "Hdf5_file::read<MArray1d>: data to read has wrong rank");
        }
        dataspace.getSimpleExtentDims(&dims[0], NULL);
        MArray1i retval(boost::extents[dims[0]]);

        DataSpace memspace(rank, &dims[0]);
        int * data_out = retval.origin();
        dataset.read(data_out, atomic_type, memspace, dataspace);
        return retval;
#endif

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

template<class Archive>
    void
    Hdf5_file::save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP(file_name)
                << BOOST_SERIALIZATION_NVP(is_open)
                << BOOST_SERIALIZATION_NVP(current_flag);
        if (is_open) {
            flush();
            copy_to_serialization_directory(file_name);
        }
    }

template<class Archive>
    void
    Hdf5_file::load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP(file_name)
                >> BOOST_SERIALIZATION_NVP(is_open)
                >> BOOST_SERIALIZATION_NVP(current_flag);
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
Hdf5_file::save<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version) const;
template
void
Hdf5_file::save<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version) const;

template
void
Hdf5_file::load<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);
template
void
Hdf5_file::load<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
