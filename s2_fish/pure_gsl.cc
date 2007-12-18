#include <iostream>

#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "array_2d.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

class GSL_random
{
    private:
    gsl_rng * rng;
    public:
    GSL_random(unsigned long int seed = 0,
        const gsl_rng_type *generator = gsl_rng_ranlxd2);
    void fill_array_unit_gaussian(double *array, int num, int stride=1);
    void fill_array_uniform(double *array, int num, int stride=1);
    ~GSL_random();
};

GSL_random::GSL_random(unsigned long int seed,
        const gsl_rng_type *generator)
{
    gsl_rng_env_setup();
    rng = gsl_rng_alloc(generator);
    gsl_rng_set(rng,seed);
}

void
GSL_random::fill_array_unit_gaussian(double *array, int num, int stride)
{
    for(int i = 0; i<num; i+=stride) {
            array[i] = gsl_ran_ugaussian_ratio_method(rng);
    }
}

void
GSL_random::fill_array_uniform(double *array, int num, int stride)
{
    for(int i = 0; i<num; i+=stride) {
            array[i] = gsl_rng_uniform(rng);
    }
}

GSL_random::~GSL_random()
{
    gsl_rng_free(rng);
}

gsl_matrix
gsl_matrix_from_Array_2d(const Array_2d<double> &array)
{
    gsl_matrix matrix;
    std::vector<int> shape = array.get_shape();
    matrix.size1= shape.at(0);
    matrix.size2 = shape.at(1);
    matrix.tda = shape.at(1);
    matrix.data = array.get_data_ptr();
    matrix.block = 0;
    matrix.owner = 0;
    
    return matrix;
}

void
array_2d_to_octave_file(const Array_2d<double> &array, const std::string filename)
{
    std::ofstream stream(filename.c_str());
    for (int i=0; i<array.get_shape()[0]; ++i) {
        for (int j=0; j<array.get_shape()[1]; ++j) {
            stream << std::setprecision(16) << array(i,j);
            if (j == array.get_shape()[1]-1) {
                stream << std::endl;
            } else {
                stream << " ";
            }
        }
    }
    stream.close();
}

int main()
{
    std::cout << "hello gsl\n";
    GSL_random gslr;
    int num_particles=10000;
    Array_2d<double> a(num_particles,6);
    //~ gslr.fill_array_unit_gaussian(a.get_data_ptr()+2,a.get_size(),6);
    gslr.fill_array_unit_gaussian(a.get_data_ptr(),a.get_size());

    std::vector<double> means(6);
    for(int j=0; j<6; ++j) means[j] = 0;
    Array_2d<double> covs(6,6);
    covs.zero_all();
    for(int n=0; n<num_particles; ++n) {
        for(int j=0; j<6; ++j) means[j] += a(n,j);
    }
    for(int j=0; j<6; ++j) means[j] *= 1.0/num_particles;
    std::cout << "means: ";
    for(int j=0; j<6; ++j) std::cout << means[j] << " ";
    std::cout << std::endl;
    for(int n=0; n<num_particles; ++n) {
        for (int i=0; i<6; ++i) {
            for (int j=0; j<=i; ++j) {
                covs(i,j) += (a(n,i)-means[i])*(a(n,j)-means[j]);
            }
        }
    }
    for (int i=0; i<6; ++i) {
        for (int j=i; j<6; ++j) {
            covs(i,j) = covs(j,i);
        }
    }
    covs.scale(1.0/num_particles);
    covs.print("covs");
    
    Array_2d<double> desired_covs(6,6);
    desired_covs.zero_all();
    for (int i=0; i<6; ++i) {
        desired_covs(i,i) = 1.0*(i+1.0);
    }
    for (int i=0; i<3; ++i) {
        desired_covs(2*i,2*i+1) = desired_covs(2*i+1,2*i) = 0.2+i/10.0;
    }
    desired_covs.print("desired_covs");

    gsl_matrix C = gsl_matrix_from_Array_2d(desired_covs);
    // C -> G
    int err =  gsl_linalg_cholesky_decomp(&C);
    gsl_matrix *G = &C;
    std::cout << "gsl_linalg_cholesky_decomp returns " << err << std::endl;
    
    gsl_matrix X = gsl_matrix_from_Array_2d(covs);
    // X -> H
    err =  gsl_linalg_cholesky_decomp(&X);
    std::cout << "gsl_linalg_cholesky_decomp returns " << err << std::endl;

    for (int i=0; i<6; ++i) {
        for (int j = i+1; j<6; ++j) {
            covs(i,j) = 0.0;
            desired_covs(i,j) = 0.0;
        }
    }
    //~ desired_covs.print("G");
    //~ covs.print("H");
    
    gsl_permutation *perm = gsl_permutation_alloc(6);
    int signum;
    err = gsl_linalg_LU_decomp(&X, perm, &signum);
    std::cout << "gsl_linalg_LU_decomp returns " << err << std::endl;
    
    gsl_matrix *Hinv = gsl_matrix_alloc(6,6);
    err = gsl_linalg_LU_invert(&X, perm, Hinv);
    std::cout << "gsl_linalg_LU_invert returns " << err << std::endl;

    gsl_matrix *A = gsl_matrix_alloc(6,6);
    
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, G, Hinv, 0.0, A);
    
    Array_2d<double> rho_sub(a);
    for(int n=0; n<num_particles; ++n) {
        for(int j=0; j<6; ++j) rho_sub(n,j) -= means[j];
    }
    
    gsl_matrix rho_sub_gsl = gsl_matrix_from_Array_2d(rho_sub);
    gsl_matrix r = gsl_matrix_from_Array_2d(a);
    
    gsl_blas_dgemm(CblasNoTrans, CblasTrans,1.0,&rho_sub_gsl,A,0.0,&r);
    // n.b.: we have not provided for a desired non-zero r, i.e., rbar.
    
    gsl_permutation_free(perm);
    gsl_matrix_free(Hinv);
    gsl_matrix_free(A);
    
    array_2d_to_octave_file(a,"r.dat");
    std::cout << "success!\n";
    return 0;
}
