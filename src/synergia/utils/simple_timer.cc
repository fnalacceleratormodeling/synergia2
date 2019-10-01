
#include "synergia/utils/simple_timer.h"


std::map<std::string, simple_timer_counter::timing> 
simple_timer_counter::timings = std::map<std::string, simple_timer_counter::timing>();


void simple_timer_print(Logger & logger)
{
#ifdef SIMPLE_TIMER
    logger(LoggerV::INFO) 
        //<< std::scientific 
        << std::setprecision(8)
        << "simple_timer:\n"
        << "(label)    \t(time)    \t(count)\n";

    for(auto const& timing : simple_timer_counter::timings)
    {
        logger 
            << timing.first << ": \t" 
            << timing.second.sum << ", \t"
            << timing.second.count << "\n";
    }

    logger << "\n";
#endif
}
