// -*- C++ -*-
#ifndef HAVE_ND_ARRAY_H
#define HAVE_ND_ARRAY_H true

#include <vector>
#include <complex>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>

template<class T>
class Nd_array;

template<class T>
struct complex_helper;

template<class T >
struct complex_helper<std::complex<T> > {
  typedef Nd_array<T> the_type;
  the_type real(Nd_array<std::complex<T> > * t);
  the_type imag(Nd_array<std::complex<T> > * t);
};

template<class T>
struct complex_helper {
  typedef Nd_array<T> the_type;
  the_type real(Nd_array<T> * t);
  the_type imag(Nd_array<T> * t);
};

template<class T>
class Nd_array {
 private:
  std::vector<T> storage;
  std::vector<int> dims;
  bool frozen;
  int vector_index(int const indices[]) const;
  void print_recursive(std::string name, int which_index, int indices[]) const;
  void set_dims(int order, int const dims_in[]);
  inline void assert_dims(int const indices[]) const;
 public:
  Nd_array();
  Nd_array(int order, int const dims_in[]);
  Nd_array(std::vector<int> const dims_in);
  Nd_array(const Nd_array& original);

  template<class T1> void copy(Nd_array<T1> * T1_tensor);
  void copy_data(T* base_address, int order, int const dims_in[]);
  // real and imag do not work (!) To be fixed during a time of great patience.
  typename complex_helper<T>::the_type real();
  typename complex_helper<T>::the_type imag();

  void reshape(int order, int const dims_in[]);
  void reshape(std::vector<int> const& dims);
  void freeze_shape();
  std::vector<int> get_shape() const;

  void zero_all( );
  void set(int const indices[],T val);
  void set(std::vector<int> const& indices,T val);
  void add_to_point(int const indices[],T val);
  void add_to_point(std::vector<int> const& indices,T val);
  T get(int const indices[]) const;
  T get(std::vector<int> const& indices) const;

  void scale(T factor);
  void add(T constant);

  int get_length();
  T* get_base_address();

  void describe() const;
  void print(std::string name) const;

  void write_to_fstream(std::ofstream& stream);
  void write_to_file(std::string filename);
  void read_from_fstream(std::ifstream& stream);
  void read_from_file(std::string filename);
};

typedef Nd_array<double> Real_nd_array;
typedef Nd_array<std::complex<double> > Complex_nd_array;

template<class T> std::vector<T> read_line_vector(std::ifstream& stream);

// Begin implementation.
#include <iostream>
#include <iomanip>

template<class T>
Nd_array<T> 
complex_helper<std::complex<T> >::real(Nd_array<std::complex<T> > * t) {
  Nd_array<T> retval(t->get_shape());
  for(int i=0; i<t->get_length(); ++i) {
  	*(retval->get_base_address()+i) = (t->get_base_address()+i)->real();
  }
  return retval;
}

template<class T>
Nd_array<T> 
complex_helper<std::complex<T> >::imag(Nd_array<std::complex<T> > * t) {
  Nd_array<T> retval(t->get_shape());
  for(int i=0; i<t->get_length(); ++i) {
  	*(retval->get_base_address()+i) = (t->get_base_address()+i)->imag();
  }
  return retval;
}

template<class T>
Nd_array<T> 
complex_helper<T>::real(Nd_array<T> * t) {
  Nd_array<T> retval(t->get_shape());
  for(int i=0; i<t->get_length(); ++i) {
  	*(retval->get_base_address()+i) = *(t->get_base_address()+i);
  }
  return retval;
}

template<class T>
Nd_array<T> 
complex_helper<T>::imag(Nd_array<T> * t) {
  Nd_array<T> retval(t->get_shape());
  for(int i=0; i<t->get_length(); ++i) {
  	*(retval->get_base_address()+i) = 0;
  }
  return retval;
}

template<class T>
Nd_array<T>::Nd_array()
{
  frozen = false;
}

template<class T>
void
Nd_array<T>::set_dims(int order, int const dims_in[])
{
  int size = 1;
  dims.resize(order);
  for (int i = 0; i < order ; i++) {
    dims[i] = dims_in[i];
    size *= dims[i];
  }
  storage.resize(size);
}

template<class T>
Nd_array<T>::Nd_array(int order, int const dims_in[])
{
  frozen = false;
  set_dims(order,dims_in);
}

template<class T>
Nd_array<T>::Nd_array(std::vector<int> const dims_in)
{
  frozen = false;
  set_dims(dims_in.size(),&dims_in[0]);
}

template<class T>
Nd_array<T>::Nd_array(const Nd_array& original)
{
  frozen = false;
  storage = original.storage;
  dims = original.dims;
}


template<class T>
template<class T1>
void 
Nd_array<T>::copy(Nd_array<T1> * T1_tensor)
{
  std::vector<int> shape(T1_tensor->get_shape());
  set_dims(shape.size(),&shape[0]);
  for(int i=0; i<T1_tensor->get_length(); ++i) {
  	storage[i] = *(T1_tensor->get_base_address()+i);
  }
}

template<class T>
void 
Nd_array<T>::copy_data(T* base_address, int order, int const dims_in[])
{
  reshape(order,dims_in);
  for(int i=0; i<get_length(); ++i) {
  	storage[i] = *(base_address+i);
  }
}

template<class T>
typename complex_helper<T>::the_type
Nd_array<T>::real()
{
  return complex_helper<T>::real(this);
}

template<class T>
typename complex_helper<T>::the_type
Nd_array<T>::imag()
{
  return complex_helper<T>::imag(this);
}

template<class T>
int
Nd_array<T>::get_length()
{
  return storage.size();
}

template<class T>
T*
Nd_array<T>::get_base_address()
{
  return &storage[0];
}

template<class T>
void
Nd_array<T>::reshape(int order, int const dims_in[])
{
  bool shape_changed = false;
  if (order == dims.size()) {
    for (int i=0; i<order; ++i) {
      if (dims[i] != dims_in[i]) {
	shape_changed = true;
      }
    }
  } else {
    shape_changed = true;
  }
  if (shape_changed && frozen) {
    throw 
      std::runtime_error("Attempt to change the shape of a frozen Nd_array");
  }
  if (shape_changed) {
    set_dims(order, dims_in);
  }
}

template<class T>
void
Nd_array<T>::reshape(std::vector<int> const& dims)
{
  reshape(dims.size(),&dims[0]);
}

template<class T>
void
Nd_array<T>::freeze_shape()
{
  frozen = true;
}

template<class T>
std::vector<int>
Nd_array<T>::get_shape() const
{
  return dims;
}


template<class T>
void
Nd_array<T>::zero_all( )
{
  int size = storage.size( );
  for ( int i = 0; i < size; ++i )
    {
      storage[i] = 0.0;
    }
}

template<class T>
int
Nd_array<T>::vector_index(int const indices[]) const
{
  int val = 0;
  int multiplier = 1;
  for(int i = 0; i < dims.size(); i++) {
    val += indices[i]*multiplier;
    multiplier *= dims[i];
  }
  return val;
}

template<class T>
T
Nd_array<T>::get(int const indices[]) const
{
  assert_dims(indices);
  return storage[vector_index(indices)];
}

template<class T>
T
Nd_array<T>::get(std::vector<int> const& indices) const
{
  assert_dims(&indices[0]);
  return storage[vector_index(&indices[0])];
}

template<class T>
void
Nd_array<T>::assert_dims(int const indices[]) const
{
  if (dims.size() == 0) {
    std::stringstream message("");
    message << "Nd_array range error: operation attempted on zero-dimensional array.";
    throw std::out_of_range(message.str());
  }
  // optimize for case size == 3
  if (dims.size() == 3) {
    if ((indices[0] < 0) || (indices[0] >= dims[0]) ||
	(indices[1] < 0) || (indices[1] >= dims[1]) ||
	(indices[2] < 0) || (indices[2] >= dims[2])) {
      std::stringstream message("");
      message << "Nd_array range error: point ("
	      << indices[0] << ","
	      << indices[1] << ","
	      << indices[2] << ") outside of bounds of field ("
	      << "0-" << dims[0]-1 << ","
	      << "0-" << dims[1]-1 << ","
	      << "0-" << dims[2]-1 << ")";
      throw std::out_of_range(message.str());
    }
  } else {
    bool out_of_bounds = false;
    for(int i=0; i<dims.size(); ++i) {
      if ((indices[i] < 0) || (indices[i] >= dims[i])) {
	out_of_bounds = true;
      }
    }
    if (out_of_bounds) {
      std::stringstream message("");
      message << "Nd_array range error: point (";
      for(int i=0; i<dims.size(); ++i) {
	message << indices[i] ;
	if (i<dims.size()-1) {
	  message << ",";
	}
      }
      message << ") outside of bounds (";
      for(int i=0; i<dims.size(); ++i) {
	message << "0-" << dims[i]-1;
	if (i<dims.size()-1) {
	  message << ",";
	}
      }
      message << ")";
      throw std::out_of_range(message.str());
    }
  }
}
 
template<class T>
void 
Nd_array<T>::set(int const indices[], T val)
{
  assert_dims(indices);
  storage[vector_index(indices)] = val;
}

template<class T>
void 
Nd_array<T>::set(std::vector<int> const& indices,T val)
{
  assert_dims(&indices[0]);
  storage[vector_index(&indices[0])] = val;
}

template<class T>
void 
Nd_array<T>::add_to_point(int const indices[], T val)
{
  assert_dims(indices);
  storage[vector_index(indices)] += val;
}

template<class T>
void 
Nd_array<T>::add_to_point(std::vector<int> const& indices,T val)
{
  assert_dims(&indices[0]);
  storage[vector_index(&indices[0])] += val;
}

template<class T>
void
Nd_array<T>::scale(T factor)
{
  int size = storage.size();
  for ( int i = 0; i < size; ++i )
    {
      storage[i] *= factor;
    }
}

template<class T>
void
Nd_array<T>::add(T constant)
{
  int size = storage.size();
  for ( int i = 0; i < size; ++i )
    {
      storage[i] += constant;
    }
}

template<class T>
void
Nd_array<T>::describe() const
{
  std::cout << "tensor has order " << dims.size() << " ";
  std::cout << "and dimensions (";
    for (int i = 0; i < dims.size(); i++) {
      if (i>0) {
	std::cout << ",";
      }
      std::cout << dims[i];
    }
  std::cout << ")\n";
}

template<class T>
void
Nd_array<T>::print_recursive(std::string name, int which_index,
			     int indices[]) const
{
  if (which_index == 2) {
    std::cout << name << "(:,:";
    for(int i=2; i<dims.size(); i++) {
      std::cout << "," << indices[i];
    }
    std::cout << ")" << std::endl;
    for(int i=0; i<dims[0]; i++) {
      indices[0] = i;
      for(int j=0; j<dims[1]; j++) {
	indices[1] = j;
	std::cout << std::setw(9);
	std::cout << get(indices);
      }
      std::cout << std::endl;
    }
  } else {
    for(int i=0; i<dims[which_index-1]; i++) {
      indices[which_index-1] = i;
      print_recursive(name,which_index-1,indices);
      std::cout << std::endl;
    }
  }
}

template<class T>
void
Nd_array<T>::print(std::string name) const
{
  if (dims.size() == 0) {
    std::cout << name << ": empty\n";
  } else if(dims.size() == 1) {
    std::cout << name << ":\n";
    int index[1];
    for (int i=0; i<dims[0]; ++i) {
      index[0] = i;
      std::cout << std::setw(9);
      std::cout << get(index) << std::endl;
    }
  } else {
    int indices[dims.size()];
    print_recursive(name,dims.size(),indices);
  }
}

template<class T>
void
Nd_array<T>::write_to_fstream(std::ofstream& stream)
{
  int size(1);
  for(int i=0; i<dims.size(); ++i) {
    size *= dims[i];
    stream << dims[i];
    if(i<dims.size()-1) {
      stream << " ";
    }
  }
  stream << std::endl;
  
  for(int i=0; i<size; ++i) {
    stream << storage[i];
    if(i<size-1) {
      stream << " ";
    }
  }
  stream << std::endl;
}

template<class T>
void
Nd_array<T>::write_to_file(std::string filename)
{
  std::ofstream stream(filename.c_str());
  write_to_fstream(stream);
  stream.close();
}

template<class T>
std::vector<T>
read_line_vector(std::vector<T>& v, std::ifstream& stream)
{
  const int buffer_length(1000);
  char buffer[buffer_length];
  stream.getline(buffer,buffer_length);
  std::stringstream sstream(buffer);
  while (! sstream.eof()) {
    T val;
    sstream >> val;
    v.push_back(val);
  }
  return v;
}

template<class T>
void
Nd_array<T>::read_from_fstream(std::ifstream& stream)
{
  std::vector<int> read_dims;
  read_line_vector(read_dims,stream);
  int size(1);
  for(int i=0; i<read_dims.size(); ++i) {
    size *= read_dims[i];
  }
  reshape(read_dims.size(),&read_dims[0]);
  
  for(int i=0; i<size; ++i) {
    stream >> storage[i];
  }
}

template<class T>
void
Nd_array<T>::read_from_file(std::string filename)
{
  std::ifstream stream(filename.c_str());
  read_from_fstream(stream);
  stream.close();
}

#endif
