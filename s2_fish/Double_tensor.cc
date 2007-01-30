#include "Double_tensor.h"
#include <iostream>
#include <iomanip>

Double_tensor::Double_tensor()
{
}

Double_tensor::Double_tensor(int order, int const dims_in[])
{
  int size = 1;
  dims.reserve(order);
  for (int i = 0; i < order ; i++) {
    dims.push_back(dims_in[i]);
    size *= dims[i];
  }
  storage.reserve(size);
}

Double_tensor::Double_tensor(const Double_tensor& original)
{
  storage = original.storage;
  dims = original.dims;
}

int
Double_tensor::length()
{
  return storage.size();
}

double*
Double_tensor::get_base_address()
{
  return &storage[0];
}

void
Double_tensor::reshape(int order, int const dims_in[])
{
  int size = 1;
  dims.resize(order);
  for (int i = 0; i < order ; i++) {
    dims[i] = dims_in[i];
    size *= dims[i];
  }
  storage.resize(size);
}

std::vector<int>
Double_tensor::get_shape()
{
  return dims;
}

int
Double_tensor::vector_index(int const indices[]) const
{
  int val = 0;
  int multiplier = 1;
  for(int i = 0; i < dims.size(); i++) {
    val += indices[i]*multiplier;
    multiplier *= dims[i];
  }
  return val;
}

double &
Double_tensor::operator() (int const indices[])
{
  return storage[vector_index(indices)];
}

double
Double_tensor::operator() (int const indices[]) const
{
  return storage[vector_index(indices)];
}

void 
Double_tensor::set(int const indices[], double val)
{
  storage[vector_index(indices)] = val;
}

std::vector<int>
Double_tensor::get_dims()
{
  return dims;
}

void
Double_tensor::describe()
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

void
Double_tensor::print_n(int which_index, int indices[]) const
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

void
Double_tensor::print() const
{
  int indices[dims.size()];
  print_n(dims.size(),indices);
}

void 
Double_tensor::zero_all( )
{
  int size = storage.size( );
  for ( int i = 0; i < size; ++i )
    {
      storage[i] = 0.0;
    }
}
