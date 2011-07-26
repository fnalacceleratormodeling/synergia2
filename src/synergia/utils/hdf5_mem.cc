#include <iostream>
#include "simple_memory_tracker.h"
#include "hdf5_serial_writer.h"

void
doit0()
{
    H5::H5File file("hdf5_mem0.h5", H5F_ACC_TRUNC);
//    Hdf5_serial_writer<int > writer(file, "i");
//    file.close();
}

void
doit1e5()
{
    H5::H5File file("hdf5_mem1e5.h5", H5F_ACC_TRUNC);
    Hdf5_serial_writer<int > writer(file, "i");
    for (int i = 0; i < 100000; ++i) {
        writer.append(i);
    }
    file.close();
}

int
main()
{
    Simple_memory_tracker memtracker;
    memtracker.print("start (absolute usage)");

    doit0();
    memtracker.print("after doit0 (relative)");

    doit0();
    memtracker.print("after doit0 again (relative)");

    doit0();
    memtracker.print("after doit0 yet again (relative)");

    doit1e5();
    memtracker.print("after doit1e5 (relative)");

    doit1e5();
    memtracker.print("after doit1e5 again (relative)");

    memtracker.print("end (relative)");
}
