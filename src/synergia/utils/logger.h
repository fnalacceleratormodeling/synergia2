#ifndef LOGGER_H_
#define LOGGER_H_
#include <ostream>
#include <fstream>
#include <string>
#include "synergia/utils/commxx.h"

class Logger
{
private:
    std::ostream * stream_ptr;
    bool have_stream;
    std::ofstream * fstream_ptr;
    bool have_fstream;
public:
    /// Log to screen on a single rank
    Logger(int rank, bool log = true);
    /// Log to file and, optionally, screen on a single rank
    Logger(int rank, std::string const& filename, bool screen = true, bool log =
            true);
    /// Log to file and, optionally, screen on a single rank
    /// This variation on the previous constructor is needed
    /// for technical C++ reasons...
    Logger(int rank, char const * filename, bool screen = true,
            bool log = true);
    /// Log to a separate file on each rank
    Logger(std::string const& filename_base, bool log = true);
    Logger &
    set_stream(std::ostream & stream);
    Logger &
    write(std::string const& str);
    template<typename T>
        Logger &
        operator<<(T const& t)
        {
            if (have_stream) {
                (*stream_ptr) << t;
            }
            if (have_fstream) {
                (*fstream_ptr) << t;
            }
            return *this;
        }
    Logger &
    operator<<(std::ostream &
    (*op)(std::ostream &));
    Logger &
    flush();
    ~Logger();
};

#endif /* LOGGER_H_ */
