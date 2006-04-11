#include <iostream>
#include <list>
#undef _POSIX_C_SOURCE
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include <Numeric/arrayobject.h>
#include <vector>
#include "bmlfactory.h"
#include "Mapping.h"
#include "Double_tensor.h"

extern "C" {
#include <sys/time.h>
}

double 
double_time()
{
  timeval t;
  gettimeofday(&t,NULL);
  return t.tv_sec + t.tv_usec/1.0e6;
}
  
using namespace boost::python;

int convert_chef_index(int impact_index)
{
  return impact_index/2+3*(impact_index%2);
}

double quickpow(double x, int i)
{
  double retval = 1.0;
  while (i--) {
    retval *= x;
  }
  return retval;
}

class Particles
{
private:
  double *data;
  int index(int row,int col) {
    return col*7 + row;
  };
public:
  Particles(numeric::array& numeric_particles) {
    data = reinterpret_cast<double *>
      (reinterpret_cast<PyArrayObject*>(numeric_particles.ptr())->data);
  };
  double & operator()(int component, int particle) {
    return data[index(component,particle)];
  };
};

class Linear_map
{
private:
  double *data;
  int index(int row, int col) {
    return row*7 + col;
  };
public:
  Linear_map(numeric::array& numeric_map) {
    data = reinterpret_cast<double *>
      (reinterpret_cast<PyArrayObject*>(numeric_map.ptr())->data);
  };
  double & operator()(int row, int col) {
    return data[index(row,col)];
  };
};

#define MANUAL_LOOP_UNROLL shockingly_yes

void apply_linear_map(numeric::array& numeric_particles, int num_particles,
		      numeric::array& numeric_map)
{
  double t0 = double_time();
  Particles particles(numeric_particles);
  Linear_map map(numeric_map);
  
  double temp[6];
  for(int part=0; part<num_particles; ++part) {
#ifdef MANUAL_LOOP_UNROLL
    temp[0] = 0.0;
    temp[0] += particles(0,part)*map(0,0);
    temp[0] += particles(1,part)*map(0,1);
    temp[0] += particles(2,part)*map(0,2);
    temp[0] += particles(3,part)*map(0,3);
    temp[0] += particles(4,part)*map(0,4);
    temp[0] += particles(5,part)*map(0,5);

    temp[1] = 0.0;
    temp[1] += particles(0,part)*map(1,0);
    temp[1] += particles(1,part)*map(1,1);
    temp[1] += particles(2,part)*map(1,2);
    temp[1] += particles(3,part)*map(1,3);
    temp[1] += particles(4,part)*map(1,4);
    temp[1] += particles(5,part)*map(1,5);

    temp[2] = 0.0;
    temp[2] += particles(0,part)*map(2,0);
    temp[2] += particles(1,part)*map(2,1);
    temp[2] += particles(2,part)*map(2,2);
    temp[2] += particles(3,part)*map(2,3);
    temp[2] += particles(4,part)*map(2,4);
    temp[2] += particles(5,part)*map(2,5);

    temp[3] = 0.0;
    temp[3] += particles(0,part)*map(3,0);
    temp[3] += particles(1,part)*map(3,1);
    temp[3] += particles(2,part)*map(3,2);
    temp[3] += particles(3,part)*map(3,3);
    temp[3] += particles(4,part)*map(3,4);
    temp[3] += particles(5,part)*map(3,5);

    temp[4] = 0.0;
    temp[4] += particles(0,part)*map(4,0);
    temp[4] += particles(1,part)*map(4,1);
    temp[4] += particles(2,part)*map(4,2);
    temp[4] += particles(3,part)*map(4,3);
    temp[4] += particles(4,part)*map(4,4);
    temp[4] += particles(5,part)*map(4,5);

    temp[5] = 0.0;
    temp[5] += particles(0,part)*map(5,0);
    temp[5] += particles(1,part)*map(5,1);
    temp[5] += particles(2,part)*map(5,2);
    temp[5] += particles(3,part)*map(5,3);
    temp[5] += particles(4,part)*map(5,4);
    temp[5] += particles(5,part)*map(5,5);

    particles(0,part) = temp[0];
    particles(1,part) = temp[1];
    particles(2,part) = temp[2];
    particles(3,part) = temp[3];
    particles(4,part) = temp[4];
    particles(5,part) = temp[5];
#else
    for(int i=0; i<6; ++i) {
      temp[i] = 0.0;
      for(int j=0; j<6; ++j) {
	temp[i] += particles(j,part)*map(i,j);
      }
    }
    for(int i=0; i<6; ++i) {
      particles(i,part) = temp[i];
    }
#endif
  }
}

class Fast_mapping_term
{
public:
  double coeff;
  int  *i;
  int order;
  Fast_mapping_term(int order) {
    i = new int[order];
    this->order = order;
  };
  Fast_mapping_term(const Fast_mapping_term& t) {
    coeff = t.coeff;
    order = t.order;
    i = new int[order]; 
    for (int j=0; j<order; ++j) {
      i[j] = t.i[j];
    }
  };
  ~Fast_mapping_term() {delete [] i;};
};

class Unit_conversion
{
private:
  double *data;
public:
  Unit_conversion(numeric::array& numeric_u) {
    data = reinterpret_cast<double *>
    (reinterpret_cast<PyArrayObject*>(numeric_u.ptr())->data);
  };
  double & operator()(int impact_index) {
    return data[impact_index];
  };
};

class Fast_mapping
{
private:
  std::vector<std::vector<std::list<Fast_mapping_term> > > terms;
  int order;
public:
  Fast_mapping(numeric::array& numeric_u, Mapping mapping);
  void apply(numeric::array& numeric_particles, int num_particles);
};

Fast_mapping::Fast_mapping(numeric::array& numeric_u, Mapping mapping) {
  Unit_conversion u(numeric_u);
  order = mapping.Weight();
  terms.resize(6);
  for (int comp_index=0; comp_index<6; ++comp_index) {
    terms.at(comp_index).resize(order);
  }      
  for (int i = 0; i<6 ; ++i) {
    int chef_i = convert_chef_index(i);
    const TJLterm<double,std::complex<double> > * jlterm;
    mapping(chef_i).resetConstIterator();
    int nterm = 0;
    while(jlterm = mapping(chef_i).stepConstIteratorPtr()){
      if ((jlterm->coefficient() == 0.0) || 
	  (jlterm->exponents().Sum() == 0)) {
	//ignore zero coefficients
      } else {
	if (jlterm->exponents().Sum() > 0) {
	  int term_order = jlterm->exponents().Sum();
	  Fast_mapping_term tmp_term(term_order);
	  tmp_term.coeff = jlterm->coefficient()*u(i);
	  int which = 0;
	  for (int index=0; index<6; ++index) {
	    int chef_index = convert_chef_index(index);
	    int expt = jlterm->exponents()(chef_index);
	    for (int count = 0; count < expt; ++count) {
	      tmp_term.i[which] = index;
	      ++which;
	    }
	    tmp_term.coeff *= 1.0/quickpow(u(index),expt);
	  }
	  terms.at(i).at(term_order-1).push_back(tmp_term);
	} else {
	  std::cerr 
	    << "Fast_mapping found something funky witha term of order "
	    <<jlterm->exponents().Sum() << std::endl;
	}
      }
    }
  } 
}

void 
Fast_mapping::apply(numeric::array& numeric_particles, int num_particles) {
  Particles particles(numeric_particles);
  double temp[6];
  for(int part=0; part<num_particles; ++part) {
    for (int i = 0; i<6; ++i) {
      temp[i] = 0.0;
      double term;
      std::list<Fast_mapping_term>::const_iterator telem;
      for (telem = terms[i][0].begin(); telem !=terms[i][0].end(); ++telem) {
	term = telem->coeff;
	term *= particles(telem->i[0],part);
	temp[i] += term;
      }
      if (order > 1) {
	for (telem = terms[i][1].begin(); telem !=terms[i][1].end(); ++telem) {
	  term = telem->coeff;
	  term *= particles(telem->i[0],part);
	  term *= particles(telem->i[1],part);
	  temp[i] += term;
	}
	if (order > 2) {
	  for (telem = terms[i][2].begin(); telem !=terms[i][2].end(); ++telem) {
	    term = telem->coeff;
	    term *= particles(telem->i[0],part);
	    term *= particles(telem->i[1],part);
	    term *= particles(telem->i[2],part);
	    temp[i] += term;
	  }
	  if (order > 3) {
	    for (int suborder = 1; suborder<order; ++suborder) {
	      for (telem = terms[i][suborder].begin();
		   telem !=terms[i][suborder].end();
		   ++telem) {
		term = telem->coeff;
		for (int j=0; j<=suborder; ++j) {
		  term *= particles(telem->i[j],part);
		}
		temp[i] += term;
	      }
	    }
	  }
	}
      }
    }
#ifdef MANUAL_LOOP_UNROLL
    particles(0,part) = temp[0];
    particles(1,part) = temp[1];
    particles(2,part) = temp[2];
    particles(3,part) = temp[3];
    particles(4,part) = temp[4];
    particles(5,part) = temp[5];
#else
    for(int i=0; i<6; ++i) {
      particles(i,part) = temp[i];
    }
#endif
  }
}
 
BOOST_PYTHON_MODULE(mappers)
{
  numeric::array::set_module_and_type("Numeric", "ArrayType");
  def("apply_linear_map",&apply_linear_map);
  class_<Fast_mapping>("Fast_mapping",init<numeric::array&, Mapping>() )
    .def("apply",&Fast_mapping::apply);
}

