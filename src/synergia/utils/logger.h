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
    Logger(int rank);
    /// Log to screen and file on a single rank
    Logger(int rank, std::string const& filename);
    /// Log to a separate file on each rank
    Logger(std::string const& filename_base);
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
