Code Conventions
================

A C++ class Foo will typically have two derived typedefs associated with it,
Foo_sptr and Foos. The former is a Boost shared_ptr to Foo. The latter is a
standard library list of Foo_sptr's. In Python, Foo and Foo_sptr may be
considered interchangeable.

1-, 2- and 3-dimensional arrays are implemented via Boost MultiArrays in C++
and Numpy arrays in Python. The Python bindings take care of the low-overhead
conversion between the two types.