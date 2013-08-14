#ifndef MULTI_ARRAY_TO_STRING_H_
#define MULTI_ARRAY_TO_STRING_H_

#include <boost/multi_array.hpp>
#include <iomanip>
#include <sstream>

template<typename T, int N>
    std::string
    _recursive_to_string(T a, std::string const & name,
            boost::array<boost::multi_array_types::index, N > &indices,
            const int which_index)
    {
        std::stringstream stream;
        if (which_index == 2) {
            stream << name << "(:,:";
            for (int i = 2; i < N; i++) {
                stream << "," << indices[i];
            }
            stream << ")" << std::endl;
            for (int i = a.index_bases()[0];
                    i < (a.index_bases()[0] + a.shape()[0]); i++) {
                indices[0] = i;
                for (int j = a.index_bases()[1];
                        j < (a.index_bases()[1] + a.shape()[1]); j++) {
                    indices[1] = j;
                    stream << std::setw(12);
                    stream << a(indices);
                }
                stream << std::endl;
            }
        } else {
            for (int i = a.index_bases()[which_index - 1];
                    i
                            < (a.index_bases()[which_index - 1]
                                    + a.shape()[which_index - 1]); i++) {
                indices[which_index - 1] = i;
                stream
                        << _recursive_to_string<T, N >(a, name, indices,
                                which_index - 1);
            }
        }
        return stream.str();
    }

template<typename T>
    std::string
    multi_array_to_string(T a, std::string const& name = "")
    {
        std::stringstream stream;
        const int N = T::dimensionality;
        boost::array<boost::multi_array_types::index, N > indices;
        if (N == 0) {
            stream << name << ": empty\n";
        } else if (N == 1) {
            stream << name << ":\n";
            for (indices[0] = 0; indices[0] < a.shape()[0]; ++indices[0]) {
                stream << std::setw(12);
                stream << a(indices) << std::endl;
            }
        } else {
            stream << _recursive_to_string<T, N >(a, name, indices, N);
        }
        return stream.str();
    }

#endif /* MULTI_ARRAY_TO_STRING_H_ */
