#include "synergia/utils/catch.hpp"

#include "synergia/utils/hdf5_file.h"
#include "synergia/utils/multi_array_typedefs.h"


TEST_CASE("hdf5_append", "[Hdf5_file_append]")
{
    CHECK(true);
}

TEST_CASE("hdf5_append_scalar", "[Hdf5_file_append]")
{
    int mpi_rank = Commxx::world_rank();
    int mpi_size = Commxx::world_size();

    std::stringstream ss;
    ss << "hdf5_append_scalar_" << mpi_size << ".h5";
    std::string fname = ss.str();

    {
        Hdf5_file file(fname, Hdf5_file::truncate, Commxx());

        // v1: [3, 4]
        file.append("v1", 3, false);
        file.append("v1", 4, false);

        // v2: [3.2, 4.2]
        file.append("v2", 3.2, false);
        file.append("v2", 4.2, false);

        // v3: [[3, 3, 3, 3], [4, 4, 4, 4]]
        file.append("v3", 3.0, true);
        file.append("v3", 4.0, true);

        // v4: [[10, 11, 12, 13], [20, 21, 22, 23]]
        file.append("v4", 10.0+mpi_rank, true);
        file.append("v4", 20.0+mpi_rank, true);
    }

    {
        Hdf5_file file(fname, Hdf5_file::read_write, Commxx());

        // v1: [3, 4]
        file.append("v1", 5, false);
        file.append("v1", 6, false);

        file.append("v3", 3.0, true);
        file.append("v3", 4.0, true);
    }

    {
        Hdf5_file file(fname, Hdf5_file::read_only, Commxx());

        // v1
        auto v1 = file.read<karray1i_row>("v1");
        REQUIRE(v1.extent(0) == 4);
        CHECK(v1(0) == 3);
        CHECK(v1(1) == 4);
        CHECK(v1(2) == 5);
        CHECK(v1(3) == 6);

        // v2
        auto v2 = file.read<karray1d_row>("v2");
        REQUIRE(v2.extent(0) == 2);
        CHECK(v2(0) == 3.2);
        CHECK(v2(1) == 4.2);
    }
}


