/*******************************************
** grid.cc
** Contains:
** 
*******************************************/
#include "scalar_field.h"

#include <stdexcept>
#include <sstream>
#include <cmath>


//--------------------------------------------------------------------
Scalar_Field::Scalar_Field(  )
{
  //int tmp_points[3] = {0, 0, 0};
  //num_points = tmp_points;
  //set_physical_params( {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0} );
}

//--------------------------------------------------------------------
Scalar_Field::Scalar_Field(int3 num_points, double3 physical_size, double3 physical_offset )
{
  set_num_points( num_points );
  set_physical_params( physical_size, physical_offset );
}

//--------------------------------------------------------------------
Scalar_Field::~Scalar_Field()
{
  ;
}

//--------------------------------------------------------------------
void Scalar_Field::set_num_points( int3 num_points )
{
  this->num_points = num_points;
  points.reshape( 3, num_points.c_array());
  update_constants( );
}

//--------------------------------------------------------------------
int3 Scalar_Field::get_num_points()
{
  std::vector<int> v(points.get_shape());
  return int3(v);
}

//--------------------------------------------------------------------
void Scalar_Field::set_physical_params( double3 physical_size, double3 physical_offset )
{
  for (int i = 0; i < 3; ++i)
    {
      this->physical_size[i] = physical_size[i];
      this->physical_offset[i] = physical_offset[i];
    }
  update_constants( );
}

//--------------------------------------------------------------------
double3 Scalar_Field::get_physical_size()
{
  return physical_size;
}

//--------------------------------------------------------------------
double3 Scalar_Field::get_physical_offset()
{
  return physical_offset;
}

//--------------------------------------------------------------------
void Scalar_Field::zero_the_points( )
{
  points.zero_all( );
}

void Scalar_Field::set_point( int3 indices, double val )
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

double Scalar_Field::get_point ( int3 indices)
{
  double retval;
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

void Scalar_Field::add_to_point( int3 indices, double val )
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

int3 Scalar_Field::get_leftmost_indices( double3 location )
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
double3 Scalar_Field::get_leftmost_offsets( double3 location )
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

void Scalar_Field::update_constants( )
{
  for ( int i = 0; i < 3; ++i )
    {
      left[i] = ( physical_offset[i] - physical_size[i]/ 2);
      h[i] = ( num_points[i] - 1.0 ) / physical_size[i];
    }
}

void Scalar_Field::print_points( )
{
  points.print();
}

int Scalar_Field::array_length( )
{
  return points.length();
}

double* Scalar_Field::array_base_address()
{
  return points.get_base_address();
}
