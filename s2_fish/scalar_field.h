// -*- C++ -*-
#ifndef HAVE_SCALAR_FIELD_H
#define HAVE_SCALAR_FIELD_H

#include "nd_array.h"
#include "triple.h"
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <cmath>
#include <vector>

template<class T>
class Scalar_field
{
private:
    Nd_array<T> points;
    Double3 physical_size;
    Double3 physical_offset;

    //derived quantities
    Double3 left;
    Double3 h;

    inline T get_delta(int location[3], int axis) const;

public:
    Scalar_field();
    Scalar_field(int num_points[3], double physical_size[3],
                 double physical_offset[3]);
    Scalar_field(int num_points[3], double physical_size[3],
                 double physical_offset[3], int dim0_lower, int dim0_upper);
    Scalar_field(std::vector<int> num_points, std::vector<double> physical_size,
                 std::vector<double> physical_offset);
    Scalar_field(std::vector<int> num_points, std::vector<double> physical_size,
                 std::vector<double> physical_offset,
                 int dim0_lower, int dim0_upper);
    ~Scalar_field();
    template<class T1> void copy(Scalar_field<T1> * T1_scalar_field);

    void set_physical_params(double physical_size[3], double physical_offset[3]);
    void set_physical_params(std::vector<double> physical_size,
                             std::vector<double> physical_offset);
    std::vector<double> get_physical_size() const;
    std::vector<double> get_physical_offset() const;
    std::vector<double> get_cell_size() const;

    inline Nd_array<T>& get_points();
    inline const Nd_array<T>& get_points() const;

    int* get_leftmost_indices(double location[3]) const;
    std::vector<int> get_leftmost_indices(std::vector<double> location) const ;
    double* get_leftmost_offsets(double location[3]) const ;
    std::vector<double> get_leftmost_offsets(std::vector<double> location) const;

    inline T get_val(double location[3]) const ;
    inline T get_val(std::vector<double> location) const ;
    inline T get_deriv(double location[3], int axis) const ;

    void write_to_fstream(std::ofstream& stream);
    void write_to_file(std::string filename);
    void read_from_fstream(std::ifstream& stream);
    void read_from_file(std::string filename);
    void describe();
};

typedef Scalar_field<double> Real_scalar_field;
typedef Scalar_field<std::complex<double> > Complex_scalar_field;

// begin implementation...

template<class T>
Scalar_field<T>::Scalar_field()
{
}

template<class T>
Scalar_field<T>::Scalar_field(int num_points[3], double physical_size[3],
                              double physical_offset[3])
{
    points.reshape(3, num_points);
    set_physical_params(physical_size, physical_offset);
}

template<class T>
Scalar_field<T>::Scalar_field(int num_points[3], double physical_size[3],
                              double physical_offset[3],
                              int dim0_lower, int dim0_upper)
{
    points.reshape(3, num_points, dim0_lower, dim0_upper);
    set_physical_params(physical_size, physical_offset);
}

template<class T>
Scalar_field<T>::Scalar_field(std::vector<int> num_points,
                              std::vector<double> physical_size,
                              std::vector<double> physical_offset)
{
    if ((num_points.size() != 3) || (physical_size.size() != 3) ||
            (physical_offset.size() != 3)) {
        throw std::invalid_argument("Arguments to Scalar_field constructors must be of length 3");
    }
    points.reshape(num_points);
    set_physical_params(physical_size, physical_offset);
}

template<class T>
Scalar_field<T>::Scalar_field(std::vector<int> num_points,
                              std::vector<double> physical_size,
                              std::vector<double> physical_offset,
                              int dim0_lower, int dim0_upper)
{
    if ((num_points.size() != 3) || (physical_size.size() != 3) ||
            (physical_offset.size() != 3)) {
        throw std::invalid_argument("Arguments to Scalar_field constructors must be of length 3");
    }
    points.reshape(num_points, dim0_lower, dim0_upper);
    set_physical_params(physical_size, physical_offset);
}

template<class T>
template<class T1>
void
Scalar_field<T>::copy(Scalar_field<T1> * T1_scalar_field)
{
    points.copy(&T1_scalar_field->points);
    set_physical_params(&T1_scalar_field->get_physical_size()[0],
                        &T1_scalar_field->get_physical_offset()[0]);
}

template<class T>
Scalar_field<T>::~Scalar_field()
{
}

template<class T>
void
Scalar_field<T>::set_physical_params(double physical_size[3],
                                     double physical_offset[3])
{
    points.freeze_shape();
    for (int i = 0; i < 3; ++i) {
        this->physical_size[i] = physical_size[i];
        this->physical_offset[i] = physical_offset[i];
        left[i] = (physical_offset[i] - physical_size[i] / 2);
        h[i] = physical_size[i] / (points.get_shape()[i] - 1.0);
    }
}

template<class T>
void
Scalar_field<T>::set_physical_params(std::vector<double> physical_size,
                                     std::vector<double> physical_offset)
{
    if ((physical_size.size() != 3) || (physical_offset.size() != 3)) {
        throw std::invalid_argument("Arguments to Scalar_field::set_physical_params must be of length 3");
    }
    set_physical_params(&physical_size[0], &physical_offset[0]);
}

template<class T>
std::vector<double>
Scalar_field<T>::get_physical_size() const 
{
    return physical_size.vector();
}

template<class T>
std::vector<double>
Scalar_field<T>::get_physical_offset() const 
{
    return physical_offset.vector();
}

template<class T>
std::vector<double>
Scalar_field<T>::get_cell_size() const 
{
    return h.vector();
}

template<class T>
inline Nd_array<T>&
Scalar_field<T>::get_points()
{
    return points;
}

template<class T>
inline const Nd_array<T>&
Scalar_field<T>::get_points() const
{
    return points;
}

template<class T>
int*
Scalar_field<T>::get_leftmost_indices(double location[3]) const
{
    Int3 leftmost_val;
    for (int i = 0; i < 3; ++i) {
        leftmost_val[i] = static_cast<int>
                          (floor((location[i] - left[i]) / h[i]));
    }
    return leftmost_val;
}

template<class T>
std::vector<int>
Scalar_field<T>::get_leftmost_indices(std::vector<double> location) const
{
    return Int3(get_leftmost_indices(&location[0])).vector();
}

template<class T>
double*
Scalar_field<T>::get_leftmost_offsets(double location[3]) const
{
    Double3 offset;
    for (int i = 0; i < 3; ++i) {
        double scaled_location =  (location[i] - left[i]) / h[i];
        offset[i] = scaled_location -
                    static_cast<int>(floor(scaled_location));
    }
    return offset;
}

template<class T>
inline T
Scalar_field<T>::get_val(double location[3]) const
{
    // Interpolate between grid points. There is no unique scheme to do this
    // in 3D, so we choose to use trilinear interpolation.
    Int3 c(get_leftmost_indices(location)); // c for corner
    Double3 f(get_leftmost_offsets(location)); // f for fractional difference
    T val = ((1.0 - f[0]) * (1.0 - f[1]) * (1.0 - f[2]) * points.get(c) +
             f[0] * (1.0 - f[1]) * (1.0 - f[2]) * points.get(Int3(c[0] + 1, c[1], c[2])) +
             (1.0 - f[0]) * f[1] * (1.0 - f[2]) * points.get(Int3(c[0], c[1] + 1, c[2])) +
             (1.0 - f[0]) * (1.0 - f[1]) * f[2] * points.get(Int3(c[0], c[1], c[2] + 1)) +
             f[0] * f[1] * (1.0 - f[2]) * points.get(Int3(c[0] + 1, c[1] + 1, c[2])) +
             f[0] * (1.0 - f[1]) * f[2] * points.get(Int3(c[0] + 1, c[1], c[2] + 1)) +
             (1.0 - f[0]) * f[1] * f[2] * points.get(Int3(c[0], c[1] + 1, c[2] + 1)) +
             f[0] * f[1] * f[2] * points.get(Int3(c[0] + 1, c[1] + 1, c[2] + 1)));
    return val;
}

template<class T>
inline T
Scalar_field<T>::get_delta(int location[3], int axis) const
{
    // return delta such that derivative is delta/h
    Int3 left(location), right(location);
    double scale = 0.5;
    --left[axis];
    if (left[axis] < 0) {
        left[axis] = 0;
        scale = 1.0;
    }
    ++right[axis];
    if (right[axis] >= points.get_shape()[axis]) {
        right[axis] = points.get_shape()[axis] - 1;
        scale = 1.0;
    }
    return scale*(points.get(right) - points.get(left));
}

template<class T>
inline T
Scalar_field<T>::get_deriv(double location[3], int axis) const
{
    // Interpolate between grid points. There is no unique scheme to do this
    // in 3D, so we choose to use trilinear interpolation.
    Int3 c(get_leftmost_indices(location)); // c for corner
    Double3 f(get_leftmost_offsets(location)); // f for fractional difference
    T val = ((1.0 - f[0]) * (1.0 - f[1]) * (1.0 - f[2]) * get_delta(c, axis) +
             f[0] * (1.0 - f[1]) * (1.0 - f[2]) * get_delta(Int3(c[0] + 1, c[1], c[2]), axis) +
             (1.0 - f[0]) * f[1] * (1.0 - f[2]) * get_delta(Int3(c[0], c[1] + 1, c[2]), axis) +
             (1.0 - f[0]) * (1.0 - f[1]) * f[2] * get_delta(Int3(c[0], c[1], c[2] + 1), axis) +
             f[0] * f[1] * (1.0 - f[2]) * get_delta(Int3(c[0] + 1, c[1] + 1, c[2]), axis) +
             f[0] * (1.0 - f[1]) * f[2] * get_delta(Int3(c[0] + 1, c[1], c[2] + 1), axis) +
             (1.0 - f[0]) * f[1] * f[2] * get_delta(Int3(c[0], c[1] + 1, c[2] + 1), axis) +
             f[0] * f[1] * f[2] * get_delta(Int3(c[0] + 1, c[1] + 1, c[2] + 1), axis)) / h[axis];
    return val;
}

template<class T>
inline T
Scalar_field<T>::get_val(std::vector<double> location) const
{
    return get_val(&location[0]);
}

template<class T>
std::vector<double>
Scalar_field<T>::get_leftmost_offsets(std::vector<double> location) const
{
    return Double3(get_leftmost_offsets(&location[0])).vector();
}

template<class T>
void
Scalar_field<T>::write_to_fstream(std::ofstream& stream)
{
    stream << physical_size[0] << " "
    << physical_size[1] << " "
    << physical_size[2] << std::endl;
    stream << physical_offset[0] << " "
    << physical_offset[1] << " "
    << physical_offset[2] << std::endl;
    points.write_to_fstream(stream);
}

template<class T>
void
Scalar_field<T>::write_to_file(std::string filename)
{
    std::ofstream stream(filename.c_str());
    write_to_fstream(stream);
    stream.close();
}

template<class T>
void
Scalar_field<T>::read_from_fstream(std::ifstream& stream)
{
    std::vector<double> read_physical_size, read_physical_offset;
    read_line_vector(read_physical_size, stream);
    read_line_vector(read_physical_offset, stream);
    points.read_from_fstream(stream);
    set_physical_params(read_physical_size, read_physical_offset);
}

template<class T>
void
Scalar_field<T>::read_from_file(std::string filename)
{
    std::ifstream stream(filename.c_str());
    read_from_fstream(stream);
    stream.close();
}

template<class T>
void
Scalar_field<T>::describe()
{
    points.describe();
}
#endif				//	ifndef HAVE_SCALAR_FIELD_H

