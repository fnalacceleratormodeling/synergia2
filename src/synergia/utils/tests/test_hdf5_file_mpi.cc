#include "synergia/utils/catch.hpp"

#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"


TEST_CASE("hdf5_file_dim_check_array_nothrow", "[Hdf5_file]")
{
    Hdf5_file file("hdf5_file_test_dim_array_nothrow.h5", 
            Hdf5_file::truncate, Commxx());

    std::vector<int> vi(6);
    std::iota(vi.begin(), vi.end(), Commxx::world_rank() * 100);

    CHECK_NOTHROW(file.write("vi", vi.data(), vi.size(), true));
    CHECK_NOTHROW(file.write("vi2", vi.data(), vi.size(), false));
}

TEST_CASE("hdf5_file_dim_check_array", "[Hdf5_file]")
{
    Hdf5_file file("hdf5_file_test_dim_array.h5", 
            Hdf5_file::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::vector<int> vi(6 + mpi_rank);
    std::iota(vi.begin(), vi.end(), mpi_rank * 100);

    CHECK_NOTHROW(file.write("vi", vi.data(), vi.size(), true));

    if (mpi_size > 1) CHECK_THROWS(file.write("vi2", vi.data(), vi.size(), false));
    else CHECK_NOTHROW(file.write("vi2", vi.data(), vi.size(), false));
}

TEST_CASE("hdf5_file_dim_check_zero_sized", "[Hdf5_file]")
{
    Hdf5_file file("hdf5_file_test_dim_array_zero_sized.h5", 
            Hdf5_file::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::vector<int> vi(mpi_rank);
    std::iota(vi.begin(), vi.end(), mpi_rank * 100);

    CHECK_NOTHROW(file.write("vi", vi.data(), vi.size(), true));

    if (mpi_size > 1) CHECK_THROWS(file.write("vi2", vi.data(), vi.size(), false));
    else CHECK_NOTHROW(file.write("vi2", vi.data(), vi.size(), false));
}


TEST_CASE("hdf5_file_dim_check_kv_nothrow", "[Hdf5_file]")
{
    Hdf5_file file("hdf5_file_test_dim_kv_nothrow.h5", 
            Hdf5_file::truncate, Commxx());

    karray2d_row arr1("arr1", 4, 3);
    arr1(2, 1) = 2.0;
    arr1(3, 2) = 3.0;

    CHECK_NOTHROW(file.write("arr1", arr1, true));
    CHECK_NOTHROW(file.write("arr2", arr1, false));
}

TEST_CASE("hdf5_file_dim_check_kv_1", "[Hdf5_file]")
{
    Hdf5_file file("hdf5_file_test_dim_kv_1.h5", 
            Hdf5_file::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    karray2d_row arr1("arr1", 4 + mpi_rank, 3);
    arr1(2, 1) = 2.0;
    arr1(3, 2) = 3.0;

    CHECK_NOTHROW(file.write("arr1", arr1, true));

    if( mpi_size > 1) CHECK_THROWS(file.write("arr2", arr1, false));
    else CHECK_NOTHROW(file.write("arr2", arr1, false));
}

TEST_CASE("hdf5_file_dim_check_kv_2", "[Hdf5_file]")
{
    Hdf5_file file("hdf5_file_test_dim_kv_2.h5", 
            Hdf5_file::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    karray2d_row arr1("arr1", 4 + mpi_rank, 3 + mpi_rank);
    arr1(2, 1) = 2.0;
    arr1(3, 2) = 3.0;

    if( mpi_size > 1) CHECK_THROWS(file.write("arr1", arr1, true));
    else CHECK_NOTHROW(file.write("arr1", arr1, true));

    if( mpi_size > 1) CHECK_THROWS(file.write("arr2", arr1, false));
    else CHECK_NOTHROW(file.write("arr2", arr1, false));
}

TEST_CASE("hdf5_file_dim_check_kv_zero_dim", "[Hdf5_file]")
{
    Hdf5_file file("hdf5_file_test_dim_kv_zero_dim.h5", 
            Hdf5_file::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    karray2d_row arr1("arr1", mpi_rank, 3);
    //arr1(2, 1) = 2.0;
    //arr1(3, 2) = 3.0;

    CHECK_NOTHROW(file.write("arr1", arr1, true));

    if( mpi_size > 1) CHECK_THROWS(file.write("arr2", arr1, false));
    else CHECK_NOTHROW(file.write("arr2", arr1, false));
}

TEST_CASE("hdf5_file_dim_check_kv_diff_rank", "[Hdf5_file]")
{
    Hdf5_file file("hdf5_file_test_dim_kv_diff_rank.h5", 
            Hdf5_file::truncate, Commxx());

    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    if (mpi_rank == 0)
    {
        karray2d_row arr1("arr1", 4, 3);

        if (mpi_size == 1)
        {
            CHECK_NOTHROW(file.write("arr1", arr1, true));
            CHECK_NOTHROW(file.write("arr2", arr1, false));
        }
        else
        {
            CHECK_THROWS(file.write("arr1", arr1, true));
            CHECK_THROWS(file.write("arr2", arr1, false));
        }
    }
    else
    {
        karray3d_row arr1("arr1", 4, 3, 2);

        CHECK_THROWS(file.write("arr1", arr1, true));
        CHECK_THROWS(file.write("arr2", arr1, false));
    }
}

TEST_CASE("hdf5_file_scalar_read", "[Hdf5_file]")
{
    {
        Hdf5_file file("hdf5_file_scalar_read.h5", 
                Hdf5_file::truncate, Commxx());

        int i = 3;
        CHECK_NOTHROW(file.write("int", i, false));
    }

    {
        Hdf5_file file("hdf5_file_scalar_read.h5", 
                Hdf5_file::read_only, Commxx());

        auto i = file.read<int>("int");

        CHECK(i == 3);
    }
}


TEST_CASE("hdf5_file_kv_read", "[Hdf5_file]")
{
    {
        Hdf5_file file("hdf5_file_kv_read.h5", 
                Hdf5_file::truncate, Commxx());

        karray2d_row arr1("arr1", 4, 3);
        arr1(0, 0) = 10.0;
        arr1(1, 0) = 11.0;
        arr1(2, 0) = 12.0;
        arr1(3, 0) = 13.0;

        arr1(1, 1) = 1.0;
        arr1(2, 1) = 2.0;
        arr1(3, 2) = 3.0;

        CHECK_NOTHROW(file.write("arr1", arr1, false));
    }

    {
        Hdf5_file file("hdf5_file_kv_read.h5", 
                Hdf5_file::read_only, Commxx());

        auto arr1 = file.read<karray2d_row>("arr1");

        REQUIRE(arr1.extent(0) == 4);
        REQUIRE(arr1.extent(1) == 3);

        CHECK(arr1(1, 1) == 1.0);
        CHECK(arr1(2, 1) == 2.0);
        CHECK(arr1(3, 2) == 3.0);
    }


    {
        int mpi_size = Commxx::world_size();
        int mpi_rank = Commxx::world_rank();

        if (mpi_size <= 4)
        {
            Hdf5_file file("hdf5_file_kv_read.h5", 
                    Hdf5_file::read_only, Commxx());

            std::vector<int> rows;

            switch(mpi_size)
            {
                case 1: rows = {4}; break;
                case 2: rows = {1, 3}; break;
                case 3: rows = {1, 1, 2}; break;
                case 4: rows = {1, 1, 0, 2}; break;
                default: break;
            }

            auto arr1 = file.read<karray2d_row>("arr1", rows[mpi_rank]);

            REQUIRE(arr1.extent(0) == rows[mpi_rank]);
            REQUIRE(arr1.extent(1) == 3);

            for(int i=0; i<rows[mpi_rank]; ++i)
            {
                int off = 0;
                for(int r=0; r<mpi_rank; ++r) off += rows[r];

                CHECK(arr1(i, 0) == (10.0 + off + i));
            }
        }
    }
}

TEST_CASE("hdf5_file", "[Hdf5_file]")
{
    {
        Hdf5_file file("hdf5_file_test.h5", Hdf5_file::truncate, Commxx());

        file.write_collective("int", 5);
        file.write_single("int2", 6);

        int mpi_rank = Commxx::world_rank();
        std::vector<int> vi(6);
        for(int i=0; i<vi.size(); ++i) vi[i] = i + mpi_rank * 100;

        file.write("vi", vi.data(), vi.size(), true);
        file.write("vi2", vi.data(), vi.size(), false);

        karray2d_row arr1("arr1", 4, 3);
        arr1(2, 1) = 2.0;
        arr1(3, 2) = 3.0;

        file.write("arr1", arr1, true);
        file.write("arr2", arr1, false);
    }

    {
        Hdf5_file file("hdf5_file_test.h5", Hdf5_file::read_only, Commxx());

        int mpi_rank = Commxx::world_rank();
        int mpi_size = Commxx::world_size();

        int sz = mpi_rank == 0 ? mpi_size*6 : 0;

        std::vector<int> vi(sz, 3);
        file.read("vi", vi.data(), vi.size());
    }
}

