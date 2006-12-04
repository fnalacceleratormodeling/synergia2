/*******************************************
** grid.h
** Contains:
** 
*******************************************/

#ifndef HAVE_SCALAR_FIELD_H
#define HAVE_SCALAR_FIELD_H

#include "Double_tensor.h"
#include <iostream>

class Scalar_Field
{
 private:
  Double_tensor points;
  double physical_size[3];
  double physical_offset[3];
  int num_points[3];

  double left[3];
  double h[3];
  
  void update_constants( );

 public:
  Scalar_Field(  );
  Scalar_Field(int num_points[3], double physical_size[3], double physical_offset[3] );
  ~Scalar_Field();

  void set_num_points( int num_points[3] );
  void set_physical_params( double physical_size[3], double physical_offset[3] );

  void zero_the_points( );
  void set_point( int indices[3], double val );
  double get_point( int indices[3] );
  void add_to_point( int indices[3], double val );
  int* get_nearest_indices( double location[3] );
};

#endif				//	ifndef GRID_H

