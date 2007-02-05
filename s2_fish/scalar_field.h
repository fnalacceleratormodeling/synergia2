/*******************************************
 ** scalar_field.h
 ** Contains:
 ** 
 *******************************************/

#ifndef HAVE_SCALAR_FIELD_H
#define HAVE_SCALAR_FIELD_H

#include "Double_tensor.h"
#include "TripleT.h"
#include <iostream>

template<class T>
class Scalar_field
{
 private:
  double3 physical_size;
  double3 physical_offset;
  int3 num_points;

  double3 left;
  double3 h;
  
  void update_constants( );

 public:
  T_tensor<T> points;
  Scalar_field(  );
  Scalar_field(int3 num_points, double3 physical_size, double3 physical_offset );
  ~Scalar_field();
  template<class T1> void copy(Scalar_field<T1> * T1_scalar_field);

  void set_num_points( int3 num_points );
  int3 get_num_points();
  void set_physical_params( double3 physical_size, double3 physical_offset );
  double3 get_physical_size();
  double3 get_physical_offset();
  double3 get_cell_size();
  T_tensor<T> get_points();
  void zero_the_points( );
  void set_point( int3 indices, T val );
  T get_point( int3 indices );
  void add_to_point( int3 indices, T val );
  int3 get_leftmost_indices( double3 location );
  double3 get_leftmost_offsets( double3 location );
  void print_points();
  int array_length();
  T* array_base_address();
};

typedef Scalar_field<double> Real_scalar_field;
typedef Scalar_field<std::complex<double> > Complex_scalar_field;
// begin implementation...
#include <stdexcept>
#include <sstream>
#include <cmath>

//--------------------------------------------------------------------
template<class T>
Scalar_field<T>::Scalar_field(  )
{
  //int tmp_points[3] = {0, 0, 0};
  //num_points = tmp_points;
  //set_physical_params( {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0} );
}

//--------------------------------------------------------------------
template<class T>
Scalar_field<T>::Scalar_field(int3 num_points, double3 physical_size, double3 physical_offset )
{
  set_num_points( num_points );
  set_physical_params( physical_size, physical_offset );
}

template<class T>
template<class T1> 
void
Scalar_field<T>::copy(Scalar_field<T1> * T1_scalar_field)
{
  set_num_points( T1_scalar_field->get_num_points() );
  set_physical_params( T1_scalar_field->get_physical_size(),
		       T1_scalar_field->get_physical_offset());
  points.copy(&T1_scalar_field->points);
}
//--------------------------------------------------------------------
template<class T>
Scalar_field<T>::~Scalar_field()
{
  ;
}

//--------------------------------------------------------------------
template<class T>
void Scalar_field<T>::set_num_points( int3 num_points )
{
  this->num_points = num_points;
  points.reshape( 3, num_points.c_array());
  update_constants( );
}

//--------------------------------------------------------------------
template<class T>
int3 Scalar_field<T>::get_num_points()
{
  std::vector<int> v(points.get_shape());
  return int3(v);
}

//--------------------------------------------------------------------
template<class T>
void Scalar_field<T>::set_physical_params( double3 physical_size, double3 physical_offset )
{
  for (int i = 0; i < 3; ++i)
    {
      this->physical_size[i] = physical_size[i];
      this->physical_offset[i] = physical_offset[i];
    }
  update_constants( );
}

//--------------------------------------------------------------------
template<class T>
double3 Scalar_field<T>::get_physical_size()
{
  return physical_size;
}

//--------------------------------------------------------------------
template<class T>
double3 Scalar_field<T>::get_physical_offset()
{
  return physical_offset;
}

//--------------------------------------------------------------------
template<class T>
double3 Scalar_field<T>::get_cell_size()
{
  std::vector<int> n(points.get_shape());
  return double3(physical_size[0]/(n[0]-1),
  		 physical_size[1]/(n[1]-1),
  		 physical_size[2]/(n[2]-1));
}

template<class T>
T_tensor<T> Scalar_field<T>::get_points()
{
  return points;
}
//--------------------------------------------------------------------
template<class T>
void Scalar_field<T>::zero_the_points( )
{
  points.zero_all( );
}

template<class T>
void Scalar_field<T>::set_point( int3 indices, T val )
{
  if ((indices[0] >= 0) && (indices[0] < num_points[0]) &&
      (indices[1] >= 0) && (indices[1] < num_points[1]) &&
      (indices[2] >= 0) && (indices[2] < num_points[2])) {
    points.set( indices.c_array(), val );
  } else {
    std::stringstream message("");
    message << "Scalar_field range error: point ("
	    << indices[0] << ","
	    << indices[1] << ","
	    << indices[2] << ") outside of bounds of field ("
	    << "0-" << num_points[0]-1 << ","
	    << "0-" << num_points[1]-1 << ","
	    << "0-" << num_points[2]-1 << ")";
    throw std::out_of_range(message.str());
  }
}

template<class T>
T Scalar_field<T>::get_point ( int3 indices)
{
  T retval;
  if ((indices[0] >= 0) && (indices[0] < num_points[0]) &&
      (indices[1] >= 0) && (indices[1] < num_points[1]) &&
      (indices[2] >= 0) && (indices[2] < num_points[2])) {
    retval = points(indices.c_array());
  }
  else {
    std::stringstream message("");
    message << "Scalar_field range error: point ("
	    << indices[0] << ","
	    << indices[1] << ","
	    << indices[2] << ") outside of bounds of field ("
	    << "0-" << num_points[0]-1 << ","
	    << "0-" << num_points[1]-1 << ","
	    << "0-" << num_points[2]-1 << ")";
    throw std::out_of_range(message.str());
  }
  return retval;
}

template<class T>
void Scalar_field<T>::add_to_point( int3 indices, T val )
{
  if ((indices[0] >= 0) && (indices[0] < num_points[0]) &&
      (indices[1] >= 0) && (indices[1] < num_points[1]) &&
      (indices[2] >= 0) && (indices[2] < num_points[2])) {
    points.set( indices.c_array(), val + points(indices.c_array()) );
  } else {
    std::stringstream message("");
    message << "Scalar_field range error: point ("
	    << indices[0] << ","
	    << indices[1] << ","
	    << indices[2] << ") outside of bounds of field ("
	    << "0-" << num_points[0]-1 << ","
	    << "0-" << num_points[1]-1 << ","
	    << "0-" << num_points[2]-1 << ")";
    throw std::out_of_range(message.str());
  }
}

template<class T>
int3 Scalar_field<T>::get_leftmost_indices( double3 location )
{
  int3 leftmost_val;
  for ( int i = 0; i < 3; ++i )
    {
      leftmost_val[i] = static_cast<int>
	( floor( h[i] * ( location[i] - left[i] ) ) );
    }
  return leftmost_val;
}

// return the fractional offset
template<class T>
double3 Scalar_field<T>::get_leftmost_offsets( double3 location )
{
  double3 offset;
  for ( int i = 0; i < 3; ++i )
    {
      double scaled_location =  h[i] * ( location[i] - left[i] );
      offset[i] = scaled_location - 
	static_cast<int>(floor(scaled_location));
    }
  return offset;
}

template<class T>
void Scalar_field<T>::update_constants( )
{
  for ( int i = 0; i < 3; ++i )
    {
      left[i] = ( physical_offset[i] - physical_size[i]/ 2);
      h[i] = ( num_points[i] - 1.0 ) / physical_size[i];
    }
}

template<class T>
void Scalar_field<T>::print_points( )
{
  points.print();
}

template<class T>
int Scalar_field<T>::array_length( )
{
  return points.length();
}

template<class T>
T* Scalar_field<T>::array_base_address()
{
  return points.get_base_address();
}

#endif				//	ifndef HAVE_SCALAR_FIELD_H

