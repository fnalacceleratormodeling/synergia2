#ifndef HAVE_DOUBLE_TENSOR_H
#define HAVE_DOUBLE_TENSOR_H true

#include <vector>
#include <complex>

template<class T>
class T_tensor;

template<class T>
struct complex_helper;

template<class T >
struct complex_helper<std::complex<T> > {
  typedef T_tensor<T> the_type;
  the_type real(T_tensor<std::complex<T> > * t);
  the_type imag(T_tensor<std::complex<T> > * t);
};

template<class T>
struct complex_helper {
  typedef T_tensor<T> the_type;
  the_type real(T_tensor<T> * t);
  the_type imag(T_tensor<T> * t);
};

template<class T>
class T_tensor {
 private:
  std::vector<T> storage;
  std::vector<int> dims;
  int vector_index(int const indices[]) const;
  void print_n(int which_index, int indices[]) const;
  void set_dims(int order, int const dims_in[]);
 public:
  T_tensor();
  T_tensor(int order, int const dims_in[]);
  T_tensor(std::vector<int> const dims_in);
  T_tensor(const T_tensor& original);
  template<class T1> void copy(T_tensor<T1> * T1_tensor);
  typename complex_helper<T>::the_type real();
  typename complex_helper<T>::the_type imag();
  int length();
  T* get_base_address();
  void reshape(int order, int const dims_in[]);
  std::vector<int> get_shape();
  T & operator() (int const indices[]);
  T operator() (int const indices[]) const;
  void set(int const indices[],T val);
  std::vector<int> get_dims();
  void describe();
  void print() const;
  void zero_all( );
};

typedef T_tensor<double> Double_tensor;
typedef T_tensor<std::complex<double> > Complex_tensor;

// Begin implementation. Stupid template mechanism....
#include <iostream>
#include <iomanip>

template<class T>
T_tensor<T> 
complex_helper<std::complex<T> >::real(T_tensor<std::complex<T> > * t) {
  T_tensor<T> retval(t->get_shape());
  for(int i=0; i<t->length(); ++i) {
  	*(retval->get_base_address()+i) = (t->get_base_address()+i)->real();
  }
  return retval;
}

template<class T>
T_tensor<T> 
complex_helper<std::complex<T> >::imag(T_tensor<std::complex<T> > * t) {
  T_tensor<T> retval(t->get_shape());
  for(int i=0; i<t->length(); ++i) {
  	*(retval->get_base_address()+i) = (t->get_base_address()+i)->imag();
  }
  return retval;
}

template<class T>
T_tensor<T> 
complex_helper<T>::real(T_tensor<T> * t) {
  T_tensor<T> retval(t->get_shape());
  for(int i=0; i<t->length(); ++i) {
  	*(retval->get_base_address()+i) = *(t->get_base_address()+i);
  }
  return retval;
}

template<class T>
T_tensor<T> 
complex_helper<T>::imag(T_tensor<T> * t) {
  T_tensor<T> retval(t->get_shape());
  for(int i=0; i<t->length(); ++i) {
  	*(retval->get_base_address()+i) = 0;
  }
  return retval;
}

template<class T>
T_tensor<T>::T_tensor()
{
}

template<class T>
void
T_tensor<T>::set_dims(int order, int const dims_in[])
{
  int size = 1;
  dims.reserve(order);
  for (int i = 0; i < order ; i++) {
    dims.push_back(dims_in[i]);
    size *= dims[i];
  }
  storage.reserve(size);
}

template<class T>
T_tensor<T>::T_tensor(int order, int const dims_in[])
{
  set_dims(order,dims_in);
}

template<class T>
T_tensor<T>::T_tensor(std::vector<int> const dims_in)
{
  set_dims(dims_in.size(),&dims_in[0]);
}

template<class T>
T_tensor<T>::T_tensor(const T_tensor& original)
{
  storage = original.storage;
  dims = original.dims;
}

template<class T>
template<class T1>
void 
T_tensor<T>::copy(T_tensor<T1> * T1_tensor)
{
  std::vector<int> shape(T1_tensor->get_shape());
  set_dims(shape.size(),&shape[0]);
  for(int i=0; i<T1_tensor->length(); ++i) {
  	*(this->get_base_address()+i) = *(T1_tensor->get_base_address()+i);
  }
}

template<class T>
typename complex_helper<T>::the_type
T_tensor<T>::real()
{
  return complex_helper<T>::real(this);
}

template<class T>
typename complex_helper<T>::the_type
T_tensor<T>::imag()
{
  return complex_helper<T>::imag(this);
}

template<class T>
int
T_tensor<T>::length()
{
  return storage.size();
}

template<class T>
T*
T_tensor<T>::get_base_address()
{
  return &storage[0];
}

template<class T>
void
T_tensor<T>::reshape(int order, int const dims_in[])
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
std::vector<int>
T_tensor<T>::get_shape()
{
  return dims;
}

template<class T>
int
T_tensor<T>::vector_index(int const indices[]) const
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
T &
T_tensor<T>::operator() (int const indices[])
{
  return storage[vector_index(indices)];
}

template<class T>
T
T_tensor<T>::operator() (int const indices[]) const
{
  return storage[vector_index(indices)];
}

template<class T>
void 
T_tensor<T>::set(int const indices[], T val)
{
  storage[vector_index(indices)] = val;
}

template<class T>
std::vector<int>
T_tensor<T>::get_dims()
{
  return dims;
}

template<class T>
void
T_tensor<T>::describe()
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
T_tensor<T>::print_n(int which_index, int indices[]) const
{
  if (which_index == 2) {
    std::cout << "(:,:";
    for(int i=2; i<dims.size(); i++) {
      std::cout << "," << indices[i];
    }
    std::cout << ")" << std::endl;
    for(int i=0; i<dims[0]; i++) {
      indices[0] = i;
      for(int j=0; j<dims[1]; j++) {
	indices[1] = j;
	std::cout << std::setw(9);
	std::cout << (*this)(indices);
      }
      std::cout << std::endl;
    }
  } else {
    for(int i=0; i<dims[which_index-1]; i++) {
      indices[which_index-1] = i;
      print_n(which_index-1,indices);
      std::cout << std::endl;
    }
  }
}

template<class T>
void
T_tensor<T>::print() const
{
  int indices[dims.size()];
  print_n(dims.size(),indices);
}

template<class T>
void
T_tensor<T>::zero_all( )
{
  int size = storage.size( );
  for ( int i = 0; i < size; ++i )
    {
      storage[i] = 0.0;
    }
}
#endif
