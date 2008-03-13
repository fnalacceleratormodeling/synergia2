#include <iostream>

#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "populate.h"

#include "math_constants.h"

gsl_rng * _saved_rng = 0;
const gsl_rng_type * _saved_rng_type = 0;

class GSL_random
{
private:
    gsl_rng * rng;
    const gsl_rng_type * rng_type;
public:
    GSL_random(bool init = true, unsigned long int seed = 0,
               const gsl_rng_type *generator = gsl_rng_ranlxd2);
    void fill_array_unit_gaussian(Array_nd<double> &array);
    void fill_array_uniform(Array_nd<double> &array, const double min,
                            const double max);
    ~GSL_random();
};

GSL_random::GSL_random(bool init, unsigned long int seed,
                       const gsl_rng_type *generator)
{
    gsl_rng_env_setup();
    rng_type = generator;
    rng = gsl_rng_alloc(rng_type);
    if (init || _saved_rng == 0 || _saved_rng_type != generator) {
        gsl_rng_set(rng, seed);
    } else {
        gsl_rng_memcpy(rng, _saved_rng);
    }
}

void
GSL_random::fill_array_unit_gaussian(Array_nd<double> &array)
{
    for (Array_nd<double>::Iterator it = array.begin();
            it != array.end();
            ++it) {
        *it = gsl_ran_ugaussian_ratio_method(rng);
    }
}

void
GSL_random::fill_array_uniform(Array_nd<double> &array,
                               const double min, const double max)
{
    for (Array_nd<double>::Iterator it = array.begin();
            it != array.end();
            ++it) {
        *it = gsl_ran_flat(rng, min, max);
    }
}

GSL_random::~GSL_random()
{
    if (_saved_rng != 0) {
        gsl_rng_free(_saved_rng);
    }
    _saved_rng = gsl_rng_alloc(rng_type);
    gsl_rng_memcpy(_saved_rng, rng);
    _saved_rng_type = rng_type;
    gsl_rng_free(rng);
}

class GSL_quasirandom
{
private:
    gsl_qrng * qrng;
    const gsl_qrng_type * qrng_type;
    int dimension;
public:
    GSL_quasirandom(const gsl_qrng_type *generator = gsl_qrng_sobol);
    void fill_array_unit_gaussian(Array_2d<double> &array);
    void fill_array_unit_gaussian_transverse(Array_2d<double> &array);
    //~ void fill_array_uniform(Array_nd<double> &array, const double min,
                            //~ const double max);
    ~GSL_quasirandom();
};

GSL_quasirandom::GSL_quasirandom(const gsl_qrng_type *generator)
{
    qrng_type = generator;
    dimension = 6;
    qrng = gsl_qrng_alloc(qrng_type,dimension);
}

void
GSL_quasirandom::fill_array_unit_gaussian(Array_2d<double> &array)
{
    std::vector<double> qrand_vec(dimension);
    for (int i=0; i<array.get_shape()[1]; ++i) {
        gsl_qrng_get(qrng,&qrand_vec[0]);
        for (int j=0; j<dimension; ++j) {
            array(j,i) =  gsl_cdf_gaussian_Pinv(qrand_vec[j],1.0);
        }
    }
}

void
GSL_quasirandom::fill_array_unit_gaussian_transverse(Array_2d<double> &array)
{
    std::vector<double> qrand_vec(dimension);
    for (int i=0; i<array.get_shape()[1]; ++i) {
        gsl_qrng_get(qrng,&qrand_vec[0]);
        for (int j=0; j<dimension; ++j) {
            if (j==4) {
                array(j,i) = qrand_vec[j];
            } else {
                array(j,i) =  gsl_cdf_gaussian_Pinv(qrand_vec[j],1.0);
            }
        }
    }
}

//~ void
//~ GSL_qrandom::fill_array_uniform(Array_nd<double> &array,
                               //~ const double min, const double max)
//~ {
    //~ for (Array_nd<double>::Iterator it = array.begin();
            //~ it != array.end();
            //~ ++it) {
        //~ *it = gsl_ran_flat(rng, min, max);
    //~ }
//~ }

GSL_quasirandom::~GSL_quasirandom()
{
    gsl_qrng_free(qrng);
}

gsl_matrix
gsl_matrix_from_Array_2d(const Array_2d<double> &array)
{
    gsl_matrix matrix;
    std::vector<int> shape = array.get_shape();
    matrix.size1 = shape.at(0);
    matrix.size2 = shape.at(1);
    matrix.tda = shape.at(1);
    matrix.data = array.get_data_ptr();
    matrix.block = 0;
    matrix.owner = 0;

    return matrix;
}

Array_2d<double>
Array_2d_from_gsl_matrix(const gsl_matrix &matrix)
{
    return Array_2d<double>(matrix.size1, matrix.size2, matrix.data);
}

void
array_2d_to_octave_file(const Array_2d<double> &array, const std::string filename)
{
    std::ofstream stream(filename.c_str());
    for (int i = 0; i < array.get_shape()[0]; ++i) {
        for (int j = 0; j < array.get_shape()[1]; ++j) {
            stream << std::setprecision(16) << array(i, j);
            if (j == array.get_shape()[1] - 1) {
                stream << std::endl;
            } else {
                stream << " ";
            }
        }
    }
    stream.close();
}

void
get_means_covariances(const Array_2d<double> &array,
                      Array_1d<double> &means,
                      Array_2d<double> &covariances)
{
    int num = array.get_shape()[1];
    means.set_all(0.0);
    covariances.set_all(0.0);
    for (int n = 0; n < num; ++n) {
        for (int j = 0; j < 6; ++j) means(j) += array(j, n);
    }
    means.scale(1.0 / num);
    for (int n = 0; n < num; ++n) {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j <= i; ++j) {
                covariances(i, j) += (array(i, n) - means(i)) * (array(j, n) - means(j));
            }
        }
    }
    for (int i = 0; i < 6; ++i) {
        for (int j = i; j < 6; ++j) {
            covariances(i, j) = covariances(j, i);
        }
    }
    covariances.scale(1.0 / num);
}

void
adjust_moments(Array_2d<double> &array_in,
               Array_2d<double> &array_out,
               const Array_1d<double> &means,
               const Array_2d<double> &covariances)
{
    // Verify that covariance matrix is symmetric
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            if (covariances.at(i, j) != covariances.at(j, i)) {
                throw
                std::runtime_error("covariance matrix in adjust_moments is not symmetric");
            }
        }
    }

    Array_1d<double> actual_means(6);
    Array_2d<double> actual_covs(6, 6);
    get_means_covariances(array_in, actual_means, actual_covs);

    // Calculate G
    gsl_matrix C = gsl_matrix_from_Array_2d(covariances);
    int err =  gsl_linalg_cholesky_decomp(&C);
    if (err != 0) {
        throw
        std::runtime_error("gsl_linalg_cholesky_decomp failed for matrix C");
    }
    gsl_matrix *G_ptr = &C;

    // Calculate H
    // actual_covs is const and cholesky overwrites, so work on a copy
    gsl_matrix X_tmp = gsl_matrix_from_Array_2d(actual_covs);
    gsl_matrix *X_ptr = gsl_matrix_alloc(6, 6);
    gsl_matrix_memcpy(X_ptr, &X_tmp);
    // X -> H
    err =  gsl_linalg_cholesky_decomp(X_ptr);
    if (err != 0) {
        throw
        std::runtime_error("gsl_linalg_cholesky_decomp failed for matrix X");
    }
    gsl_matrix *H_ptr = X_ptr;

    // H_ptr and G_ptr are currently G + transp(G) and H + transp(H), so
    // we need to zero out the transp() parts
    for (int i = 0; i < 6; ++i) {
        for (int j = i + 1; j < 6; ++j) {
            gsl_matrix_set(H_ptr, i, j, 0.0);
            gsl_matrix_set(G_ptr, i, j, 0.0);
        }
    }

    // Calculate H^-1
    gsl_permutation *perm = gsl_permutation_alloc(6);
    int signum;
    err = gsl_linalg_LU_decomp(H_ptr, perm, &signum);
    if (err != 0) {
        throw
        std::runtime_error("gsl_linalg_LU_decomp  failed for matrix H");
    }
    gsl_matrix *Hinv_ptr = gsl_matrix_alloc(6, 6);
    err = gsl_linalg_LU_invert(H_ptr, perm, Hinv_ptr);
    if (err != 0) {
        throw
        std::runtime_error("gsl_linalg_LU_invert  failed for matrix H");
    }

    // Calculate A
    gsl_matrix *A = gsl_matrix_alloc(6, 6);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, G_ptr, Hinv_ptr, 0.0, A);

    // Subtract actual means
    for (int n = 0; n < array_in.get_shape()[1]; ++n) {
        for (int j = 0; j < 6; ++j) array_in(j, n) -= actual_means(j);
    }

    gsl_matrix array_in_gsl = gsl_matrix_from_Array_2d(array_in);
    gsl_matrix array_out_gsl = gsl_matrix_from_Array_2d(array_out);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, A, &array_in_gsl, 0.0, &array_out_gsl);

    gsl_matrix_free(X_ptr);
    gsl_permutation_free(perm);
    gsl_matrix_free(Hinv_ptr);
    gsl_matrix_free(A);
}

void
populate_6d_gaussian(Array_2d<double> &particles,
                     const Array_1d<double> &means, const Array_2d<double> &covariances,
                     const int id_offset, const unsigned long int seed, bool init_generator)
{
    if (particles.get_shape()[0] != 7) {
        throw
        std::runtime_error("populate_6d_gaussian expects a particle array with shape (num_particles,7)");
    }
    int num_particles = particles.get_shape()[1];
    // Use the memory allocated for particles as a scratch area until the very end
    Array_2d<double> tmp(6, num_particles, particles.get_data_ptr());
    Array_2d<double> tmp2(6, num_particles);

    GSL_random gslr(init_generator, seed);
    gslr.fill_array_unit_gaussian(tmp);

    adjust_moments(tmp, tmp2, means, covariances);

    // Fill output array, adding (requested) means and setting particle ID.
    for (int n = 0; n < num_particles; ++n) {
        for (int j = 0; j < 6; ++j) particles(j, n) = tmp2(j, n) + means(j);
        particles(6, n) = (n + id_offset) * 1.0;
    }
}

void
populate_6d_gaussian_quasi(Array_2d<double> &particles,
                     const Array_1d<double> &means, const Array_2d<double> &covariances,
                     const int id_offset)
{
    if (particles.get_shape()[0] != 7) {
        throw
        std::runtime_error("populate_6d_gaussian_quasi expects a particle array with shape (num_particles,7)");
    }
    int num_particles = particles.get_shape()[1];
    // Use the memory allocated for particles as a scratch area until the very end
    Array_2d<double> tmp(6, num_particles, particles.get_data_ptr());
    Array_2d<double> tmp2(6, num_particles);

    GSL_quasirandom gslqr;
    gslqr.fill_array_unit_gaussian(tmp);

    adjust_moments(tmp, tmp2, means, covariances);

    // Fill output array, adding (requested) means and setting particle ID.
    for (int n = 0; n < num_particles; ++n) {
        for (int j = 0; j < 6; ++j) particles(j, n) = tmp2(j, n) + means(j);
        particles(6, n) = (n + id_offset) * 1.0;
    }
}

void
populate_transverse_gaussian(Array_2d<double> &particles,
                             const Array_1d<double> &means, const Array_2d<double> &covariances,
                             const int id_offset, const unsigned long int seed, bool init_generator)
{
    if (particles.get_shape()[0] != 7) {
        throw
        std::runtime_error("populate_6d_gaussian expects a particle array with shape (num_particles,7)");
    }
    int num_particles = particles.get_shape()[1];
    // Use the memory allocated for particles as a scratch area until the very end
    Array_2d<double> tmp(6, num_particles, particles.get_data_ptr());
    Array_2d<double> tmp2(6, num_particles);

    // It is simplest to let z be gaussian now, then replace it later
    GSL_random gslr(init_generator, seed);
    gslr.fill_array_unit_gaussian(tmp);

    // Symmetry requires no correlations with the z coordinate. Make a copy
    // of the covariance matrix and manually set all correlations to zero.
    Array_2d<double> covariances_modified(covariances);
    covariances_modified.copy();
    for (int k = 0; k < 6; ++k) {
        covariances_modified(4, k) = covariances_modified(k, 4) = 0.0;
    }
    covariances_modified(4, 4) = 1.0;

    adjust_moments(tmp, tmp2, means, covariances_modified);

    // Fill output array, adding (requested) means and setting particle ID.
    for (int n = 0; n < num_particles; ++n) {
        for (int j = 0; j < 6; ++j) particles(j, n) = tmp2(j, n) + means(j);
        particles(6, n) = (n + id_offset) * 1.0;
    }
    // Real z distribution
    Array_1d<double> z = particles.slice(vector2(Range(4), Range()));
    gslr.fill_array_uniform(z, -pi, pi);
}

void
populate_transverse_gaussian_quasi(Array_2d<double> &particles,
                             const Array_1d<double> &means, const Array_2d<double> &covariances,
                             const int id_offset)
{
    if (particles.get_shape()[0] != 7) {
        throw
        std::runtime_error("populate_6d_gaussian_quasi expects a particle array with shape (num_particles,7)");
    }
    int num_particles = particles.get_shape()[1];
    // Use the memory allocated for particles as a scratch area until the very end
    Array_2d<double> tmp(6, num_particles, particles.get_data_ptr());
    Array_2d<double> tmp2(6, num_particles);

    GSL_quasirandom gslqr;
    gslqr.fill_array_unit_gaussian_transverse(tmp);

    // Symmetry requires no correlations with the z coordinate. Make a copy
    // of the covariance matrix and manually set all correlations to zero.
    Array_2d<double> covariances_modified(covariances);
    covariances_modified.copy();
    for (int k = 0; k < 6; ++k) {
        covariances_modified(4, k) = covariances_modified(k, 4) = 0.0;
    }
    covariances_modified(4, 4) = 1.0;

    adjust_moments(tmp, tmp2, means, covariances_modified);
    double min_z = 1.0e30;
    double max_z = -1.0e30;
    for (int n=0; n < num_particles; ++n) {
        if (tmp2(4,n) > max_z) {
            max_z = tmp2(4,n);
        }
        if (tmp2(4,n) < min_z) {
            min_z = tmp2(4,n);
        }
    }
    // Fill output array, adding (requested) means and setting particle ID.
    for (int n = 0; n < num_particles; ++n) {
        for (int j = 0; j < 6; ++j) {
            if (j == 4) {
                particles(j, n) = (tmp2(j, n)-min_z)/(max_z - min_z)*2.0*pi - pi;
            } else {
                particles(j, n) = tmp2(j, n) + means(j);
            }
        }
        particles(6, n) = (n + id_offset) * 1.0;
    }
}
