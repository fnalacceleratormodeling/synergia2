#include "s2_diagnostics.h"
#include "mpi.h"

class Vector
{
public:
    double *data;
    Vector(numeric::array& numeric_vector)
    {
        data = reinterpret_cast<double *>
               (reinterpret_cast<PyArrayObject*>(numeric_vector.ptr())->data);
    };
    inline double & operator()(int index)
    {
        return data[index];
    };
};

class Matrix
{
private:
    double *data;
    int cols;
    inline int index(int row, int col)
    {
        return row*cols + col;
    };
public:
    Matrix(numeric::array& numeric_matrix, int cols)
    {
        data = reinterpret_cast<double *>
               (reinterpret_cast<PyArrayObject*>(numeric_matrix.ptr())->data);
        this->cols = cols;
    };
    inline double & operator()(int row, int col)
    {
        return data[index(row,col)];
    };
};

void
get_spatial_means_stds(Macro_bunch_store& mbs,
                       numeric::array& numeric_means,
                       numeric::array& numeric_stds)
{
    Vector means(numeric_means);
    Vector stds(numeric_stds);
    double sum1[3];
    double sum2[3];
    for (int i = 0; i< 3; ++i) {
        sum1[i] = 0.0;
        sum2[i] = 0.0;
    }
    for (int n = 0; n < mbs.local_num; ++n) {
        sum1[0] += mbs.local_particles(0, n);
        sum1[1] += mbs.local_particles(2, n);
        sum1[2] += mbs.local_particles(4, n);
    }
    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Allreduce(sum1,means.data,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    } else {
        for (int i=0; i<3; ++i) {
            means(i) = sum1[i];
        }
    }
    for (int i=0; i<3; ++i) {
        means(i) /= mbs.total_num;
    }
    double s[3];
    for (int n = 0; n < mbs.local_num; ++n) {
        for (int i=0; i<3; ++i) {
            s[i] = mbs.local_particles(2*i,n) - means(i);
            sum2[i] += s[i]*s[i];
        }
    }
    if (size > 1) {
        MPI_Allreduce(sum2,stds.data,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    } else {
        for (int i=0; i<3; ++i) {
            stds(i) = sum2[i];
        }
    }
    for (int i=0; i<3; ++i) {
        stds(i) = sqrt(stds(i)/mbs.total_num);
    }
}

void
get_moments_corrs(Macro_bunch_store& mbs,
                    numeric::array& numeric_units,
                    numeric::array& numeric_means,
                    numeric::array& numeric_mom2s,
                    numeric::array& numeric_corrs,
                    numeric::array& numeric_diagmom4s)
{
    Vector units(numeric_units);
    Vector means(numeric_means);
    double sum1[6];
    for (int i = 0; i< 6; ++i) {
        sum1[i] = 0.0;
    }
    for (int n = 0; n < mbs.local_num; ++n) {
        for (int i=0; i<6; ++i) {
            sum1[i] += mbs.local_particles(i, n);
        }
    }
    int rank,size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1) {
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Allreduce(sum1,means.data,6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    } else {
        for (int i=0; i<6; ++i) {
            means(i) = sum1[i];
        }
    }
    for (int i=0; i<6; ++i) {
        means(i) /= mbs.total_num;
    }
    Matrix mom2s(numeric_mom2s,6);
    Matrix corrs(numeric_corrs,6);
    Vector diagmom4s(numeric_diagmom4s);
    const int buffer_len = 36+6; // sum of lengths of mom2s and diagmom4s
    double local_buffer[buffer_len], global_buffer[buffer_len];
    for (int i=0; i<buffer_len; ++i) {
        local_buffer[i] = 0.0;
    }
    double *sum2 = local_buffer;
    double *sum4 = local_buffer + 36;
    double s[6];
    for (int n = 0; n < mbs.local_num; ++n) {
        for (int i=0; i<6; ++i) {
            s[i] = mbs.local_particles(i,n) - means(i);
            sum4[i] += s[i]*s[i]*s[i]*s[i];
            for (int j=0; j<=i; ++j) {
                sum2[6*i + j] += s[i]*s[j];
            }
        }
    }
    if (size > 1) {
        MPI_Allreduce(local_buffer,global_buffer,buffer_len,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    } else {
        for (int i=0; i<buffer_len; ++i) {
            global_buffer[i] = local_buffer[i];
        }
    }
    for (int i=0; i<6; ++i) {
        means(i) /= units(i);
    }
    for (int i=0; i<6; ++i) {
        for(int j=0; j<=i; ++j) {
            //~ mom2s(i,j) = mom2s(j,i) = global_buffer[6*i + j]/(units(i)*units(j)*mbs.total_num);
            mom2s(i,j) = global_buffer[6*i + j]/(units(i)*units(j)*mbs.total_num);
            mom2s(j,i) = mom2s(i,j);
        }
        diagmom4s(i) = global_buffer[36+i]/(units(i)*units(i)*mbs.total_num);
    }
    for (int i=0; i<6; ++i) {
        for(int j=i; j<6; ++j) {
            corrs(i,j) = corrs(j,i) = mom2s(i,j)/sqrt(mom2s(i,i)*mom2s(j,j));
        }
    }
}
