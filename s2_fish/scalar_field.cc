/*******************************************
** grid.cc
** Contains:
** 
*******************************************/
#include "scalar_field.h"

//--------------------------------------------------------------------
//	Allocate the arrays for the grid coordinates.

Scalar_Field::Scalar_Field(  )
{
  ;
}

//--------------------------------------------------------------------
//	Free Scalar_Field arrays.

Scalar_Field::~Scalar_Field()
{
	for (int j=0; j<=128; j++)
	{
	  //delete[] bc1[j];
	}
	//delete[] bc1;
}


