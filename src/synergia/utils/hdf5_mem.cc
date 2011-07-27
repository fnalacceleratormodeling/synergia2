#include <iostream>
#include "simple_memory_tracker.h"
#include "hdf5_serial_writer.h"

void
doit0(const char *filename)
{
    H5::H5File file("hdf5_mem0.h5", H5F_ACC_TRUNC);
    Hdf5_serial_writer<int > writer(file, "i");
    file.close();
}

void
doit1e5(const char *filename)
{
    H5::H5File file(filename, H5F_ACC_TRUNC);
    Hdf5_serial_writer<int > writer(file, "i");
    for (int i = 0; i < 100000; ++i) {
        writer.append(i);
    }
    file.close();
}

int
real_main(Simple_memory_tracker & memtracker)
{

    doit0("10.h5");
    memtracker.print("after doit0 (relative)");

    doit0("11.h5");
    memtracker.print("after doit0 again (relative)");

    doit0("12.h5");
    memtracker.print("after doit0 yet again (relative)");

    doit1e5("20.h5");
    memtracker.print("after doit1e5 (relative)");

    doit1e5("21.h5");
    memtracker.print("after doit1e5 again (relative)");

    doit1e5("22.h5");
    memtracker.print("after doit1e5 again (relative)");

    doit1e5("23.h5");
    memtracker.print("after doit1e5 again (relative)");


    return 0;
}

int
main()
{
    Simple_memory_tracker memtracker;
    memtracker.print("start (absolute usage)");
    int retval = real_main(memtracker);
    memtracker.print("end (absolute)", false);
    return retval;
}
