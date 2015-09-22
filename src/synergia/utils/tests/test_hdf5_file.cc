#define BOOST_TEST_MAIN
#include <algorithm>
#include <boost/test/unit_test.hpp>
#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/multi_array_check_equal.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct)
{
    Hdf5_file file("hdf5_file_empty.h5", Hdf5_file::truncate);
}

BOOST_AUTO_TEST_CASE(write_data)
{
    Hdf5_file file("hdf5_file_write.h5", Hdf5_file::truncate);
    int int_data = 7;
    file.write(int_data, "int_data");
    double double_data = 2.71828;
    file.write(double_data, "double_data");
    const int dim = 6;
    MArray1d a1d(boost::extents[dim]);
    for (int j = 0; j < dim; ++j) {
        a1d[j] = 1.1 * j;
    }
    file.write(a1d, "a1d");
    const int dim1 = 2;
    const int dim2 = 3;
    const int dim3 = 4;
    MArray3d a3d(boost::extents[dim1][dim2][dim3]);
    for (int j = 0; j < dim1; ++j) {
        for (int k = 0; k < dim2; ++k) {
            for (int l = 0; l < dim2; ++l) {
                a3d[j][k][l] = 1.1 * l + 10 * k + 100 * j;
            }
        }
    }
    file.write(a3d, "a3d");
}

namespace
{
    const char * filename = "hdf5_file.h5";
    const char * int_label = "int_data";
    const char * double_label = "double_data";
    const char * array1d_label = "array1d_data";
    const char * array2d_label = "array2d_data";
    const char * array2dfo_label = "array2dfo_data";
    const char * array3d_label = "array3d_data";
    const int num_members = 5;
}

struct Hdf5_file_fixture
{
    Hdf5_file_fixture() :
        int_data(7),
        double_data(2.71828),
        dim1(2),
        dim2(3),
        dim3(4),
        a1d(boost::extents[dim1]),
        a2d(boost::extents[dim1][dim2]),
        a3d(boost::extents[dim1][dim2][dim3])
    {
        for (int j = 0; j < dim1; ++j) {
            a1d[j] = 100 * j;
            for (int k = 0; k < dim2; ++k) {
                a2d[j][k] = 10 * k + 100 * j;
                for (int l = 0; l < dim2; ++l) {
                    a3d[j][k][l] = 1.1 * l + 10 * k + 100 * j;
                }
            }
        }
        Hdf5_file write_file(filename, Hdf5_file::truncate);
        write_file.write(int_data, int_label);
        write_file.write(double_data, double_label);
        write_file.write(a1d, array1d_label);
        write_file.write(a2d, array2d_label);
        write_file.write(a2dfo, array2dfo_label);
        write_file.write(a3d, array3d_label);
        write_file.write(a3dfo, array3dfo_label);
    }
    ~Hdf5_file_fixture()
    {
<<<<<<< HEAD
        Hdf5_file read_file(filename, Hdf5_file::read_only);
        int int_read = read_file.read<int > (int_label);
        BOOST_CHECK_EQUAL(int_read, int_data);
        double double_read = read_file.read<double > (double_label);
        BOOST_CHECK_CLOSE(double_read, double_data, tolerance);
        MArray1d a1d_read(read_file.read<MArray1d > (array1d_label));
        multi_array_check_equal(a1d_read, a1d, tolerance);
        MArray2d a2d_read(read_file.read<MArray2d > (array2d_label));
        BOOST_CHECK(a2d_read.storage_order() == a2d.storage_order());
        multi_array_check_equal(a2d_read, a2d, tolerance);
        MArray2d a2dfo_read(read_file.read<MArray2d > (array2dfo_label));
        BOOST_CHECK(a2dfo_read.storage_order() == a2dfo.storage_order());
        multi_array_check_equal(a2dfo_read, a2dfo, tolerance);
        MArray3d a3d_read(read_file.read<MArray3d > (array3d_label));
        BOOST_CHECK(a3d_read.storage_order() == a3d.storage_order());
        multi_array_check_equal(a3d_read, a3d, tolerance);
        MArray3d a3dfo_read(read_file.read<MArray3d > (array3dfo_label));
        BOOST_CHECK(a3dfo_read.storage_order() == a3dfo.storage_order());
        multi_array_check_equal(a3dfo_read, a3dfo, tolerance);
=======
>>>>>>> compacc/devel-bgq
    }

    int int_data;
    double double_data;
    int dim1;
    int dim2;
    int dim3;
    MArray1d a1d;
    MArray2d a2d;
    MArray3d a3d;
};


const double tolerance = 1.0e-13;
BOOST_FIXTURE_TEST_CASE(read_data, Hdf5_file_fixture)
{
    Hdf5_file read_file(filename, Hdf5_file::read_only);
    int int_read = read_file.read<int > (int_label);
    BOOST_CHECK_EQUAL(int_read, int_data);
    double double_read = read_file.read<double > (double_label);
    BOOST_CHECK_CLOSE(double_read, double_data, tolerance);
    MArray1d a1d_read(read_file.read<MArray1d > (array1d_label));
    multi_array_check_equal(a1d_read, a1d, tolerance);
    MArray2d a2d_read(read_file.read<MArray2d > (array2d_label));
    multi_array_check_equal(a2d_read, a2d, tolerance);
    MArray3d a3d_read(read_file.read<MArray3d > (array3d_label));
    multi_array_check_equal(a3d_read, a3d, tolerance);
}

BOOST_FIXTURE_TEST_CASE(get_member_names, Hdf5_file_fixture)
{
    Hdf5_file read_file(filename, Hdf5_file::read_only);
    std::vector<std::string > member_names(read_file.get_member_names());
    BOOST_CHECK_EQUAL(member_names.size(), num_members);
    BOOST_CHECK_EQUAL(std::count(member_names.begin(),
                                 member_names.end(), int_label), 1);
    BOOST_CHECK_EQUAL(std::count(member_names.begin(),
                                 member_names.end(), double_label), 1);
    BOOST_CHECK_EQUAL(std::count(member_names.begin(),
                                 member_names.end(), array1d_label), 1);
    BOOST_CHECK_EQUAL(std::count(member_names.begin(),
                                 member_names.end(), array2d_label), 1);
    BOOST_CHECK_EQUAL(std::count(member_names.begin(),
                                 member_names.end(), array3d_label), 1);
}

BOOST_FIXTURE_TEST_CASE(get_atomic_type, Hdf5_file_fixture)
{
    Hdf5_file read_file(filename, Hdf5_file::read_only);
    BOOST_CHECK_EQUAL(read_file.get_atomic_type(int_label), Hdf5_file::int_type);
    BOOST_CHECK_EQUAL(read_file.get_atomic_type(double_label), Hdf5_file::double_type);
    BOOST_CHECK_EQUAL(read_file.get_atomic_type(array1d_label), Hdf5_file::double_type);
    BOOST_CHECK_EQUAL(read_file.get_atomic_type(array2d_label), Hdf5_file::double_type);
    BOOST_CHECK_EQUAL(read_file.get_atomic_type(array3d_label), Hdf5_file::double_type);
}

BOOST_FIXTURE_TEST_CASE(get_dims, Hdf5_file_fixture)
{
    Hdf5_file read_file(filename, Hdf5_file::read_only);

    {
        std::vector<int > dims(read_file.get_dims(int_label));
        BOOST_CHECK_EQUAL(dims.size(), 0);
    }

    {
        std::vector<int> dims(read_file.get_dims(double_label));
        BOOST_CHECK_EQUAL(dims.size(), 0);
    }

    {
        std::vector<int> dims(read_file.get_dims(array1d_label));
        BOOST_CHECK_EQUAL(dims.size(), 1);
        BOOST_CHECK_EQUAL(dims.at(0), dim1);
    }

    {
        std::vector<int> dims(read_file.get_dims(array2d_label));
        BOOST_CHECK_EQUAL(dims.size(), 2);
        BOOST_CHECK_EQUAL(dims.at(0), dim1);
        BOOST_CHECK_EQUAL(dims.at(1), dim2);
    }

    {
        std::vector<int> dims(read_file.get_dims(array3d_label));
        BOOST_CHECK_EQUAL(dims.size(), 3);
        BOOST_CHECK_EQUAL(dims.at(0), dim1);
        BOOST_CHECK_EQUAL(dims.at(1), dim2);
        BOOST_CHECK_EQUAL(dims.at(2), dim3);
    }
}

BOOST_AUTO_TEST_CASE(test_serialize)
{
    const std::string file_name("hdf5_file_serialized.h5");
    const std::string serialize_file_name("hdf5_file.xml");
    {
        Hdf5_file file(file_name, Hdf5_file::truncate);
        file.write(1, "one");
        xml_save<Hdf5_file > (file, serialize_file_name.c_str());
    }

    {
        Hdf5_file file_resumed;
        xml_load<Hdf5_file > (file_resumed, serialize_file_name.c_str());
        BOOST_CHECK_EQUAL(1,file_resumed.read<int>("one"));
        file_resumed.write(2, "two");
    }

    Hdf5_file final_file(file_name, Hdf5_file::read_only);
    BOOST_CHECK_EQUAL(1,final_file.read<int>("one"));
    BOOST_CHECK_EQUAL(2,final_file.read<int>("two"));
}
