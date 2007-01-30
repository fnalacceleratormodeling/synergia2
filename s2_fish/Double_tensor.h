#ifndef HAVE_DOUBLE_TENSOR_H
#define HAVE_DOUBLE_TENSOR_H true

#include <vector>

class Double_tensor {
 private:
  std::vector<double> storage;
  std::vector<int> dims;
  int vector_index(int const indices[]) const;
  void print_n(int which_index, int indices[]) const;
 public:
  Double_tensor();
  Double_tensor(int order, int const dims_in[]);
  Double_tensor(const Double_tensor& original);
  int length();
  double* get_base_address();
  void reshape(int order, int const dims_in[]);
  std::vector<int> get_shape();
  double & operator() (int const indices[]);
  double operator() (int const indices[]) const;
  void set(int const indices[],double val);
  std::vector<int> get_dims();
  void describe();
  void print() const;
  void zero_all( );
};

  
#endif
