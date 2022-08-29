#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/bunch/diagnostics.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

#include "bunch_sptr_fixture.h"

const double tolerance = 1.0e-11;

BOOST_FIXTURE_TEST_CASE(construct, Bunch_sptr_fixture)
{
    Diagnostics_particles diagnostics("dummy.h5");
}


BOOST_FIXTURE_TEST_CASE(write_, Bunch_sptr_fixture)
{
    Diagnostics_particles diagnostics("diagnostics_particles_mpi1.h5");
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();

    sleep(1);

    if (comm_sptr->get_rank() == 0)
    {
        int ranks = comm_sptr->get_size();

        Hdf5_file file("diagnostics_particles_mpi1_0000.h5", Hdf5_file::read_only);
        MArray2d parts = file.read<MArray2d>("particles");

        BOOST_CHECK_EQUAL( parts.shape()[0], total_num );
        BOOST_CHECK_EQUAL( parts.shape()[1], 7 );

        for (int r=0; r<ranks; ++r)
        {
            for (int p=0; p<local_num; ++p)
            {
                for (int i=0; i<6; ++i)
                {
                    double v = 10 * p + (1.0 + p*p/1000.0) * i + 1000.0 * r;
                    BOOST_CHECK_CLOSE(parts[r*local_num + p][i], v, tolerance);
                }
            }
        }

        file.close();
    }
}

BOOST_FIXTURE_TEST_CASE(write_several, Bunch_sptr_fixture)
{
    Diagnostics_particles diagnostics("diagnostics_particles_mpi2.h5");
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();

    diagnostics.update();
    diagnostics.write();

    diagnostics.update();
    diagnostics.write();

    diagnostics.update();
    diagnostics.write();

    sleep(1);

    // check first file
    if (comm_sptr->get_rank() == 0)
    {
        int ranks = comm_sptr->get_size();

        Hdf5_file file = Hdf5_file("diagnostics_particles_mpi2_0000.h5", Hdf5_file::read_only);
        MArray2d parts = file.read<MArray2d>("particles");

        BOOST_CHECK_EQUAL( parts.shape()[0], total_num );
        BOOST_CHECK_EQUAL( parts.shape()[1], 7 );

        for (int r=0; r<ranks; ++r)
        {
            for (int p=0; p<local_num; ++p)
            {
                for (int i=0; i<6; ++i)
                {
                    double v = 10 * p + (1.0 + p*p/1000.0) * i + 1000.0 * r;
                    BOOST_CHECK_CLOSE(parts[r*local_num + p][i], v, tolerance);
                }
            }
        }

        file.close();
    }

    // check second file
    if (comm_sptr->get_rank() == 0)
    {
        int ranks = comm_sptr->get_size();

        Hdf5_file file = Hdf5_file("diagnostics_particles_mpi2_0001.h5", Hdf5_file::read_only);
        MArray2d parts = file.read<MArray2d>("particles");

        BOOST_CHECK_EQUAL( parts.shape()[0], total_num );
        BOOST_CHECK_EQUAL( parts.shape()[1], 7 );

        for (int r=0; r<ranks; ++r)
        {
            for (int p=0; p<local_num; ++p)
            {
                for (int i=0; i<6; ++i)
                {
                    double v = 10 * p + (1.0 + p*p/1000.0) * i + 1000.0 * r;
                    BOOST_CHECK_CLOSE(parts[r*local_num + p][i], v, tolerance);
                }
            }
        }

        file.close();
    }

    // check last file
    if (comm_sptr->get_rank() == 0)
    {
        int ranks = comm_sptr->get_size();

        Hdf5_file file = Hdf5_file("diagnostics_particles_mpi2_0003.h5", Hdf5_file::read_only);
        MArray2d parts = file.read<MArray2d>("particles");

        BOOST_CHECK_EQUAL( parts.shape()[0], total_num );
        BOOST_CHECK_EQUAL( parts.shape()[1], 7 );

        for (int r=0; r<ranks; ++r)
        {
            for (int p=0; p<local_num; ++p)
            {
                for (int i=0; i<6; ++i)
                {
                    double v = 10 * p + (1.0 + p*p/1000.0) * i + 1000.0 * r;
                    BOOST_CHECK_CLOSE(parts[r*local_num + p][i], v, tolerance);
                }
            }
        }

        file.close();
    }
}

BOOST_FIXTURE_TEST_CASE(write_min_max, Bunch_sptr_fixture)
{
    Diagnostics_particles diagnostics(
            "diagnostics_particles_mpi_35.h5", 3, 5);
    diagnostics.set_bunch_sptr(bunch_sptr);
    diagnostics.update();
    diagnostics.write();

    sleep(1);

    // check last file
    if (comm_sptr->get_rank() == 0)
    {
        int ranks = comm_sptr->get_size();

        Hdf5_file file = Hdf5_file("diagnostics_particles_mpi_35_0000.h5", Hdf5_file::read_only);
        MArray2d parts = file.read<MArray2d>("particles");

        BOOST_CHECK_EQUAL( parts.shape()[0], 3 );
        BOOST_CHECK_EQUAL( parts.shape()[1], 7 );

        int r = 0;

        for (int p=3; p<6; ++p)
        {
            for (int i=0; i<6; ++i)
            {
                double v = 10 * p + (1.0 + p*p/1000.0) * i + 1000.0 * r;
                BOOST_CHECK_CLOSE(parts[p-3][i], v, tolerance);
            }
        }

        file.close();
    }

}
