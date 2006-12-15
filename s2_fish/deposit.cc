#include "deposit.h"
#include <boost/python/numeric.hpp>
#include <boost/python/tuple.hpp>

// Deposit charge using Cloud-in-Cell (CIC) algorithm.
double
deposit_charge_cic(Scalar_Field& sf, Macro_bunch_store& mbs)
{
  int3 indices;
  double3 offsets;
  double total_charge = 0.0;
  for(int n=0; n<mbs.local_num; ++n) {
    double3 location(mbs.local_particles(0,n),
		     mbs.local_particles(2,n),
		     mbs.local_particles(4,n));
    indices = sf.get_leftmost_indices(location);
    offsets = sf.get_leftmost_offsets(location);
    for(int i=0; i<2; ++i) {
      for(int j=0; j<2; ++j) {
	for(int k=0; k<2; ++k) {
	  double weight = (1-i - (1-2*i)*offsets[0])* 
	    (1-j - (1-2*j)*offsets[1])* 
	    (1-k - (1-2*k)*offsets[2]);
	  try {
	    sf.add_to_point(int3(indices[0]+i,indices[1]+j,indices[2]+k),
			    weight);
	    total_charge += weight;
	  } catch(std::out_of_range e) {
	  }
	}
      }
    }
  }
  return total_charge;
}

// Deposit charge on the nearest grid point to the particle
double
deposit_charge_ngp(Scalar_Field& sf, Macro_bunch_store& mbs)
{
  int3 indices;
  double3 offsets;
  double total_charge = 0.0;
  double weight = 1.0;
  for(int n=0; n<mbs.local_num; ++n) {
    double3 location(mbs.local_particles(0,n),
		     mbs.local_particles(2,n),
		     mbs.local_particles(4,n));
    indices = sf.get_leftmost_indices(location);
    offsets = sf.get_leftmost_offsets(location);
    for(int i=0; i<3; ++i) {
      if (offsets[i]>0.5) {
	indices[i] += 1;
      }
    }
    try {
      sf.add_to_point(indices,weight);
      total_charge += weight;
    } catch(std::out_of_range e) {
    }
  }
  return total_charge;
}

