#define BOOST_TEST_MAIN
#include <vector>
#include <cstdlib>
#include <boost/test/unit_test.hpp>
#include "synergia/foundation/diagnostics_write_helper.h"
#include "synergia/utils/hdf5_serial_writer.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)

BOOST_AUTO_TEST_CASE(construct)
{
    std::string filename("test_write_helper.h5");
    bool serial = true;
    int write_skip = 1;
    Commxx commxx;
    Diagnostics_write_helper d(filename, serial, write_skip, commxx);
}

BOOST_AUTO_TEST_CASE(two_pieces)
{
    std::string filename("test_write_helper2.h5");
    bool serial = true;
    int write_skip = 1;
    Commxx commxx;
        Diagnostics_write_helper d(filename, serial, write_skip, commxx);

        Hdf5_serial_writer<double > writer_s(d.get_hdf5_file_sptr(), "s");
        for (int i = 0; i < 5; ++i) {
            double val = i * 0.01;
            writer_s.append(val);
            d.finish_write();
        }
        d.get_hdf5_file_sptr()->close();

//        Diagnostics_write_helper d2(filename, serial, write_skip, commxx);
//        d2.reopen_file();
//        Hdf5_serial_writer<double > writer_s2(d2.get_file(), "s", true);
//        for (int i = 5; i < 10; ++i) {
//            double val = i * 0.01;
//            writer_s2.append(val);
//            d2.finish_write();
//        }
//    }
}

#if 0
BOOST_AUTO_TEST_CASE(one_piece)
{
    std::string filename("test_write_helper3.h5");
    bool serial = true;
    int write_skip = 1;
    Commxx commxx;
    Diagnostics_write_helper d(filename, serial, write_skip, commxx);

    Hdf5_serial_writer<double > writer_s(d.get_file(), "s");
    for (int i = 0; i < 5; ++i) {
        double val = i * 0.01;
        writer_s.append(val);
        d.finish_write();
    }
    writer_s.close_file();
    d.close_file();
    std::string command("cp ");
    std::cout << "jfa: start copy\n";
    command += filename + " copy_" + filename;
    std::cout << "jfa: finish copy\n";
    std::system(command.c_str());
    d.reopen_file();
    writer_s.update_file(d.get_file());
    //    Hdf5_serial_writer<double > writer_s2(d.get_file(), "s", true);
    for (int i = 5; i < 10; ++i) {
        double val = i * 0.01;
        writer_s.append(val);
        d.finish_write();
    }
}
#endif
