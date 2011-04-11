#include "diagnostics_writer.h"
#include <sstream>
#include <iomanip>

//void
//Diagnostics_writer::open_file_and_init()
//{
//    std::stringstream sstream;
//    sstream << filename_base;
//    if (!diagnostics_sptr->is_serial()) {
//        sstream << "_";
//        sstream << std::setw(4);
//        sstream << std::setfill('0');
//        sstream << count;
//    }
//    sstream << filename_suffix;
//    file = H5Fcreate(sstream.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
//            H5P_DEFAULT);
//    diagnostics_sptr->init_writers(file);
//
//}
//Diagnostics_writer::Diagnostics_writer(std::string const& filename,
//        Diagnostics_sptr diagnostics_sptr) :
//    diagnostics_sptr(diagnostics_sptr), dummy(false), count(0)
//{
//    int idx = filename.rfind('.');
//    if (idx == std::string::npos) {
//        filename_base = filename;
//        filename_suffix = "";
//    } else {
//        filename_base = filename.substr(0, idx);
//        filename_suffix = filename.substr(idx);
//    }
//    if (diagnostics_sptr->is_serial()) {
//        open_file_and_init();
//    }
//}
//
//Diagnostics_writer::Diagnostics_writer() :
//    dummy(true), count(0)
//{
//}
//
//bool
//Diagnostics_writer::is_dummy() const
//{
//    return dummy;
//}
//
//Diagnostics_sptr
//Diagnostics_writer::get_diagnostics_sptr()
//{
//    return diagnostics_sptr;
//}
//
//int
//Diagnostics_writer::get_count() const
//{
//    return count;
//}
//
//void
//Diagnostics_writer::set_count(int count)
//{
//    this->count = count;
//}
//
//void
//Diagnostics_writer::write()
//{
//    if (!dummy) {
//        if (!diagnostics_sptr->is_serial()) {
//            open_file_and_init();
//        }
//        diagnostics_sptr->write_hdf5();
//        if (!diagnostics_sptr->is_serial()) {
//            H5Fclose(file);
//        }
//    }
//    ++count;
//}
//
//void
//Diagnostics_writer::update_and_write(Bunch const& bunch)
//{
//    if (!dummy) {
//        if (!diagnostics_sptr->is_serial()) {
//            open_file_and_init();
//        }
//        diagnostics_sptr->update(bunch);
//        diagnostics_sptr->write_hdf5();
//        if (!diagnostics_sptr->is_serial()) {
//            H5Fclose(file);
//        }
//    }
//    ++count;
//}
//
//Diagnostics_writer::~Diagnostics_writer()
//{
//    if (!dummy) {
//        if (diagnostics_sptr->is_serial()) {
//            H5Fclose(file);
//        }
//    }
//}
//

Multi_diagnostics::Multi_diagnostics() :
    writers()
{
}

void
Multi_diagnostics::append(
        Diagnostics_sptr diagnostics_sptr)
{
    writers.push_back(diagnostics_sptr);
}

void
Multi_diagnostics::push_back(
        Diagnostics_sptr diagnostics_sptr)
{
    writers.push_back(diagnostics_sptr);
}

Multi_diagnostics::iterator
Multi_diagnostics::begin()
{
    return writers.begin();
}

Multi_diagnostics::iterator
Multi_diagnostics::end()
{
    return writers.end();
}

Multi_diagnostics
no_diagnostics()
{
    return Multi_diagnostics();
}
