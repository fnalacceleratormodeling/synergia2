#ifndef CONTAINERS_TO_STRING_H_
#define CONTAINERS_TO_STRING_H_

#include <sstream>

template<typename T>
    std::string
    container_to_string(T const& c, const char * seperator = ",")
    {
        std::ostringstream retval;
        retval << "[";
        for (typename T::const_iterator it = c.begin(); it != c.end();) {
            retval << *it;
            ++it;
            if (it != c.end()) {
                retval << seperator;
            }
        }
        retval << "]";
        return retval.str();
    }

template<typename Tkey, typename Tval>
    std::string
    map_to_string(std::map<Tkey, Tval> const& c, const char * seperator = ",")
    {
        std::ostringstream retval;
        retval << "{";
        for (typename std::map<Tkey, Tval>::const_iterator it = c.begin(); it != c.end();) {
            retval << it->first << ":" << it->second;
            ++it;
            if (it != c.end()) {
                retval << seperator;
            }
        }
        retval << "}";
        return retval.str();
    }

template<typename Tkey>
    std::string
    map_to_string(std::map<Tkey, std::vector<double> > const& c, const char * seperator = ",")
    {
        std::ostringstream retval;
        retval << "{";
        for (typename std::map<Tkey, std::vector<double> >::const_iterator it = c.begin(); it != c.end();) {
            retval << it->first << ":" << container_to_string(it->second);
            ++it;
            if (it != c.end()) {
                retval << seperator;
            }
        }
        retval << "}";
        return retval.str();
    }


#endif /* CONTAINERS_TO_STRING_H_ */
