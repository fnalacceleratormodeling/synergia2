/*******************************************
** grid.cc
** Contains:
** 
*******************************************/
#include "scalar_field.h"

int round( double x )
{
  return static_cast<int>(x + 0.5);
}

//--------------------------------------------------------------------
Scalar_Field::Scalar_Field(  )
{
  //int tmp_points[3] = {0, 0, 0};
  //num_points = tmp_points;
  //set_physical_params( {1.0, 1.0, 1.0}, {0.0, 0.0, 0.0} );
}

//--------------------------------------------------------------------
Scalar_Field::Scalar_Field(int num_points[3], double physical_size[3], double physical_offset[3] )
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
void Scalar_Field::set_num_points( int num_points[3] )
{
  points.reshape( 3, num_points);
  update_constants( );
}

//--------------------------------------------------------------------
void Scalar_Field::set_physical_params( double physical_size[3], double physical_offset[3] )
{
  for (int i = 0; i < 3; ++i)
    {
      this->physical_size[i] = physical_size[i];
      this->physical_offset[i] = physical_offset[i];
    }
  update_constants( );
}

//--------------------------------------------------------------------
void Scalar_Field::zero_the_points( )
{
  points.zero_all( );
}

void Scalar_Field::set_point( int indices[3], double val )
{
  points.set( indices, val );
}

double Scalar_Field::get_point ( int indices[3])
{
  return points(indices);
}

void Scalar_Field::add_to_point( int indices[3], double val )
{
  points.set( indices, val + points(indices) );
}

int* Scalar_Field::get_nearest_indices( double location[3] )
{
  int nearest_val[3];
  for ( int i = 0; i < 3; ++i )
    {
      nearest_val[i] = round( h[i] * ( location[i] - left[i] ) );
    }
}

void Scalar_Field::update_constants( )
{
  for ( int i = 0; i < 3; ++i )
    {
      left[i] = ( physical_offset[i] - physical_size[i] ) / 2;
      h[i] = ( num_points[i] - 1.0 ) / physical_size[i];
    }
}
