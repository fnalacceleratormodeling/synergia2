/*
 * standalone allreduce_benchmark
 * doesn't depend on anything
 *
 * build with:
 * mpicc -o sa_allreduce_benchmark sa_allreduce_benchmark.c -lm
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include <mpi.h>

int iterations = 100000;
int gridnum = 32;

int main(int argc, char *argv[])
{
    int nranks;
    int myrank;
    int totaln;
    double *array;
    int i, j, k;
    int iter;
    double t0, t1;
    double tsum = 0.0;
    double tsum2 = 0.0;
    int parse_args(int ac, char *av[]);
    double tmean, tstd;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    printf("hello world from rank %d/%d\n", myrank, nranks);
    fflush(stdout);

    if (parse_options(argc-1, &argv[1]) != 0) return -20;

    if (myrank == 0) {
        printf("iterations: %d\n", iterations);
        printf("size: %d cubed\n", gridnum);
    }

    totaln = gridnum*gridnum*gridnum;
    if ((array = (double *) calloc(totaln, sizeof(double))) == NULL) {
        fprintf(stderr, "error rank %d: can't calloc array of %d doubles\n", myrank, totaln);
        return -10;
    }

    for (i = 0; i < gridnum; ++i) {
        for (j = 0; j < gridnum; ++j) {
            for (k = 0; k < gridnum; ++k) {
                array[i*gridnum*gridnum + j*gridnum + k] = myrank * 1.0e9 + i*1.0e6 + j*1.0e3 + k;
            }
        }
    }

    for (iter=0; iter<iterations; ++iter) {
        t0 = MPI_Wtime();
        MPI_Allreduce(MPI_IN_PLACE, (void *) array,
                      totaln, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        tsum += t1-t0;
        tsum2 += (t1-t0)*(t1-t0);
        if (myrank == 0) {
            printf("simple_timer:allreduce:%-16.9e\n", t1-t0);
            fflush(stdout);
        }
    }

    MPI_Finalize();

    tmean = tsum/iterations;
    tstd = sqrt(tsum2 - tmean*tmean*iterations);
    if (myrank == 0) {
        printf("simple_timer:allreduce.sum():%-16.9e\n", tsum);
        printf("simple_timer:allreduce.sumstd():%-16.9e\n", tstd);
        /*printf("simple_timer:allreduce.mean():%-16.9e\n", tmean); */
        /* printf("simple_timer:allreduce.std():%-16.9e\n", tstd); */
    }

    return 0;
}

int parse_options(int ac, char *av[])
{
    char *cp;

    while (ac > 0) {
        if (strstr(*av, "-h") ||
            strstr(*av, "?") ||
            strstr(*av, "help")) {
            printf("usage: sa_allreduce_benchmark [iterations=n] [gridnum=m]\n");
            return -1;
        }

        if ((cp = index(*av, '=')) == NULL) {
            fprintf(stderr, "error: unrecognized argument, should have =: %s\n", *av);
            return -1;
        }
        *cp = '\0';
        if (strcmp(*av, "iterations") == 0) {
            iterations = atoi(cp+1);
        } else if (strcmp(*av, "gridnum") == 0) {
            gridnum = atoi(cp+1);
        } else {
            fprintf(stderr, "error: invalid options: %s\n", *av);
            return -1;
        }
        --ac;
        ++av;
    }
    return 0;
}
