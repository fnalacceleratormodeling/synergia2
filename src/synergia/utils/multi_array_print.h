#ifndef MULTI_ARRAY_PRINT_H_
#define MULTI_ARRAY_PRINT_H_

#include <boost/multi_array.hpp>
#include <iostream>
#include <iomanip>
#include <string>

template<typename T, int N>
    void
    _recursive_print(T a, std::string const & name, boost::array<
            boost::multi_array_types::index, N > &indices,
            const int which_index)
    {
        if (which_index == 2) {
            std::cout << name << "(:,:";
            for (int i = 2; i < N; i++) {
                std::cout << "," << indices[i];
            }
            std::cout << ")" << std::endl;
            for (int i = a.index_bases()[0]; i < (a.index_bases()[0]
                    + a.shape()[0]); i++) {
                indices[0] = i;
                for (int j = a.index_bases()[1]; j < (a.index_bases()[1]
                        + a.shape()[1]); j++) {
                    indices[1] = j;
                    if (indices[1] > a.index_bases()[1]) {
                        std::cout << " ";
                    }
                    std::cout << std::setw(12);
                    std::cout << a(indices);
                }
                std::cout << std::endl;
            }
        } else {
            for (int i = a.index_bases()[which_index - 1]; i
                    < (a.index_bases()[which_index - 1] + a.shape()[which_index
                            - 1]); i++) {
                indices[which_index - 1] = i;
                _recursive_print<T, N > (a, name, indices, which_index - 1);
            }
        }

    }

template<typename T>
    void
    multi_array_print(T a, std::string const& name)
    {
        const int N = T::dimensionality;
        boost::array<boost::multi_array_types::index, N > indices;
        if (N == 0) {
            std::cout << name << ": empty\n";
        } else if (N == 1) {
            std::cout << name << ":\n";
            for (indices[0] = 0; indices[0] < a.shape()[0]; ++indices[0]) {
                std::cout << std::setw(12);
                std::cout << a(indices) << std::endl;
            }
        } else {
            _recursive_print<T, N > (a, name, indices, N);
        }
    }

template<typename T>
    void
    vector_print(std::vector<T > const& v, std::string const& name)
    {
        std::cout << "(";
        for (int i = 0; i < v.size(); ++i) {
            std::cout << v[i];
            if (i < (v.size() - 1)) {
                std::cout << ",";
            }
        }
        std::cout << ")";
    }

#endif /* MULTI_ARRAY_PRINT_H_ */
