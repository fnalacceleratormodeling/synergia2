#ifndef MULTI_ARRAY_ASSERT_H_
#define MULTI_ARRAY_ASSERT_H_

#include <stdexcept>
#include <sstream>

inline
void
multi_array_assert_size(Const_MArray1d_ref & a, MArray1d::size_type size,
        const char * failure_message_prefix)
{
    if (a.shape()[0] != size) {
        std::stringstream message;
        message << failure_message_prefix;
        message << " must have length ";
        message << size;
        message << ". Found length ";
        message << a.shape()[0];
        message << ".";

        throw std::runtime_error(message.str());
    }
}

inline
void
multi_array_assert_size(Const_MArray2d_ref & a, MArray2d::size_type size0,
        MArray2d::size_type size1, const char * failure_message_prefix)
{
    if ((a.shape()[0] != size0) || (a.shape()[1] != size1)) {
        std::stringstream message;
        message << failure_message_prefix;
        message << " must have size ";
        message << "(" << size0 << ", " << size1 << ")";
        message << ". Found size ";
        message << "(" << a.shape()[0] << ", " << a.shape()[1] << ")";
        message << ".";

        throw std::runtime_error(message.str());
    }
}

#endif /* MULTI_ARRAY_ASSERT_H_ */
