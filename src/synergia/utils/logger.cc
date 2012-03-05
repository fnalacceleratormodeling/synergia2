#include "logger.h"
#include <sstream>
#include "digits.h"

Logger::Logger(int rank) :
    stream_ptr(&std::cout), have_stream(false), fstream_ptr(0), have_fstream(false)
{
    if (Commxx().get_rank() == rank) {
        have_stream = true;
    }
}

Logger::Logger(int rank, std::string const& filename) :
    stream_ptr(&std::cout), have_stream(false), fstream_ptr(0), have_fstream(false)
{
    if (Commxx().get_rank() == rank) {
        have_stream = true;
        fstream_ptr = new std::ofstream(filename.c_str());
        have_fstream = true;
    }
}

Logger::Logger(std::string const& filename_base) :
    stream_ptr(&std::cout), have_stream(false), fstream_ptr(0), have_fstream(false)
{
    Commxx commxx;
    std::stringstream sstream;
    sstream << filename_base;
    sstream << "_";
    sstream << std::setw(digits(commxx.get_size()));
    sstream << std::setfill('0');
    sstream << commxx.get_rank();
    fstream_ptr = new std::ofstream(sstream.str().c_str());
    have_fstream = true;
}

Logger &
Logger::set_stream(std::ostream & stream)
{
    stream_ptr = &stream;
    return *this;
}

Logger &
Logger::write(std::string const& str)
{
    return (*this) << str;
}

Logger &
Logger::operator<<(std::ostream&
(*op)(std::ostream&))
{
    if (have_stream) {
        (*op)(*stream_ptr);
    }
    if (have_fstream) {
        (*op)(*fstream_ptr);
    }
    return *this;
}

Logger &
Logger::flush()
{
    if (have_stream) {
        stream_ptr->flush();
    }
    if (have_fstream) {
        fstream_ptr->flush();
    }
    return *this;
}

Logger::~Logger()
{
    if (have_fstream) {
        fstream_ptr->close();
        delete fstream_ptr;
    }
}
