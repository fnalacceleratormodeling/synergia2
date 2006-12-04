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
  
 public:
  Scalar_Field(  );
  ~Scalar_Field();

  
};


#endif				//	ifndef GRID_H

