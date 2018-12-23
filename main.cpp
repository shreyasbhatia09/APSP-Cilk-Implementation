//@Copyright: Jesmin Jahan Tithi, Rezaul Chowdhury, Department of Computer
// Science, Stony Brook University, Ny-11790
// Contact: jtithi@cs.stonybrook.edu,

// compile with :icpc -O3 -AVX  -xhost -o  fwr FW_rec.cpp  -ansi-alias -ip
// -opt-subscript-in-range -opt-prefetch=4 -fomit-frame-pointer
// -funroll-all-loops -vec-report  -parallel -restrict

/*
    Copy: scp -r tealab@130.245.165.107:/home/tealab/Zafar/stampede/MPI/SPAA_impl/* .
    Compile: mpicxx -O3 -AVX  -xhost -o FW_r-way FW_r-way.cpp -I$TACC_PAPI_INC -Wl,-rpath,$TACC_PAPI_LIB -L$TACC_PAPI_LIB -lpapi
    Run: ibrun -n 2 exec 8192 64
         mpiexec -n 16 exec 8192 64
*/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/timeb.h>
#include <iostream>
#include <pthread.h>
#include <math.h>
#include <immintrin.h>
#include <cstring>
#include <algorithm> 
#include <limits> // std::numeric_limits
#include <cmath>

//#include "cilktime.h"
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <papi.h>
#include <mpi.h>

using namespace std;
#ifndef TYPE
#define TYPE float
#endif

#ifndef ALIGNMENT
#define ALIGNMENT 64
#endif

// MPI_INT
#ifndef MPI_TYPE
#define MPI_TYPE MPI_FLOAT
#endif

#define ROW 0
#define COL 1
#define INF 100007
#define ERR 1e-5

// #define _DEBUG
// #define _VERIFY

void FW_APSP(TYPE *X, int n_local, int r, MPI_Comm comm_2d,
             MPI_Comm comm_row, MPI_Comm comm_col, int my_2d_rank, int my_coords[2]);

void LOG(TYPE * my_coords, char * msg);
void copy(TYPE *input, TYPE *output_1, int n);

int N, B, NP;
TYPE *dist, *X, *D;

#define min(a, b) (a < b ? a : b)
#define max(a, b) (a > b ? a : b)

void conv_RM_2_ZM_RM(TYPE *x, int ix, int jx, int ilen, int jlen) {
    if (ilen <= 0 || jlen <= 0)
        return;
    if (ilen <= B && jlen <= B) {
        for (int i = ix; i < ix + ilen; i++) {
#pragma ivdep
            for (int j = jx; j < jx + jlen; j++) {
                (*x++) = dist[(i) * (N) + j];
            }
        }
    } else {
        int n = (ilen > jlen) ? ilen : jlen;
        int c = 1;
        while (c < n)
            c = (c << 1);
        register int nn = c >> 1;

        int n11, n12, n21;
        int ni, nii, nj, njj;

        ni = min(nn, ilen);
        nj = min(nn, jlen);
        nii = max(0, ilen - nn);
        njj = max(0, jlen - nn);

        n11 = ni * nj;
        n12 = ni * njj;
        n21 = nii * nj;

        TYPE *x12, *x21, *x22;
        cilk_spawn conv_RM_2_ZM_RM(x, ix, jx, ni, nj);

        x12 = x + n11;
        cilk_spawn conv_RM_2_ZM_RM(x12, ix, jx + nj, ni, njj);

        x21 = x12 + (n12);

        cilk_spawn conv_RM_2_ZM_RM(x21, ix + ni, jx, nii, nj);

        x22 = x21 + (n21);

        conv_RM_2_ZM_RM(x22, ix + ni, jx + nj, nii, njj);

        cilk_sync;
    }
}
void conv_ZM_2_RM_RM(TYPE *x, TYPE *V, int ix, int jx, int ilen, int jlen) {
    if (ilen <= 0 || jlen <= 0)
        return;
    if (ilen <= B && jlen <= B) {
        for (int i = ix; i < ix + ilen; i++) {
#pragma ivdep
            for (int j = jx; j < jx + jlen; j++) {
                V[(i) * (N) + j] = (*x++);
            }
        }
    } else {
        int n = (ilen > jlen) ? ilen : jlen;
        int c = 1;
        while (c < n)
            c = (c << 1);
        register int nn = c >> 1;

        int n11, n12, n21;
        int ni, nii, nj, njj;

        ni = min(nn, ilen);
        nj = min(nn, jlen);
        nii = max(0, ilen - nn);
        njj = max(0, jlen - nn);

        n11 = ni * nj;
        n12 = ni * njj;
        n21 = nii * nj;

        TYPE *x12, *x21, *x22;
        cilk_spawn conv_ZM_2_RM_RM(x, V, ix, jx, ni, nj);

        x12 = x + n11;
        cilk_spawn conv_ZM_2_RM_RM(x12, V, ix, jx + nj, ni, njj);

        x21 = x12 + (n12);

        cilk_spawn conv_ZM_2_RM_RM(x21, V, ix + ni, jx, nii, nj);

        x22 = x21 + (n21);

        conv_ZM_2_RM_RM(x22, V, ix + ni, jx + nj, nii, njj);

        cilk_sync;
    }
}

/* Only for power of 2*/

void func_D(TYPE *x, TYPE *u, TYPE *v, int n) {
    if (n <= B) {
        __declspec(align(64)) TYPE V[n * n];

        int uin = 0;

#pragma parallel
        for (int i = 0; i < n; i++) {
            int in = i * n;
#pragma ivdep
            for (int j = 0; j < n; j++) {
                V[j * n + i] = v[in + j];
            }
        }

        TYPE *uu, *vv;
        int in = 0;
        uin = 0;
#pragma parallel
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                TYPE x_ij = x[in + j];

                uu = u + uin;
                vv = V + j * n;
#pragma ivdep
                for (int k = 0; k < n; k++) {
                    x_ij = min(x_ij, (*uu + *vv));
                    uu++;
                    vv++;
                }
                x[in + j] = x_ij;
            }
            in = in + n;
            uin = uin + n;
        }
        return;

    } else {
        int nn = (n >> 1);
        int nn2 = nn * nn;

        const int m11 = 0;
        int m12 = m11 + nn2;
        int m21 = m12 + nn2;
        int m22 = m21 + nn2;

        cilk_spawn func_D(x + m11, u + m11, v + m11, nn);
        cilk_spawn func_D(x + m12, u + m11, v + m12, nn);
        cilk_spawn func_D(x + m21, u + m21, v + m11, nn);
        func_D(x + m22, u + m21, v + m12, nn);

        cilk_sync;

        cilk_spawn func_D(x + m11, u + m12, v + m21, nn);
        cilk_spawn func_D(x + m12, u + m12, v + m22, nn);
        cilk_spawn func_D(x + m21, u + m22, v + m21, nn);
        func_D(x + m22, u + m22, v + m22, nn);

        cilk_sync;
    }
}

void func_C(TYPE *x, TYPE *u, TYPE *v, int n) {
    if (n <= B) {
        register TYPE uv;
        TYPE *xx, *uu, *vv;

        xx = x;
        uu = u;
        vv = v;

        for (int k = 0; k < n; k++) {
#pragma parallel
            for (int i = 0; i < n; i++) {
#pragma ivdep
                for (int j = 0; j < n; j++) {
                    uv = (*uu) + (*vv++);
                    *xx = min(uv, *xx);
                    xx++;
                }
                vv = vv - B;
                uu = uu + B;
            }
            xx = x;
            vv = vv + B;
            uu = u + k + 1;
        }
        return;
    } else {
        int nn = (n >> 1);
        int nn2 = nn * nn;

        const int m11 = 0;
        int m12 = m11 + nn2;
        int m21 = m12 + nn2;
        int m22 = m21 + nn2;

        cilk_spawn func_C(x + m11, u + m11, v + m11, nn);
        func_C(x + m21, u + m21, v + m11, nn);
        cilk_sync;
        cilk_spawn func_D(x + m12, u + m11, v + m12, nn);
        func_D(x + m22, u + m21, v + m12, nn);

        cilk_sync;

        cilk_spawn func_C(x + m12, u + m12, v + m22, nn);
        func_C(x + m22, u + m22, v + m22, nn);
        cilk_sync;

        cilk_spawn func_D(x + m11, u + m12, v + m21, nn);
        func_D(x + m21, u + m22, v + m21, nn);

        cilk_sync;
    }
}

void func_B(TYPE *x, TYPE *u, TYPE *v, int n) {
    if (n <= B) {
        register TYPE uv;
        TYPE *xx, *uu, *vv;

        xx = x;
        uu = u;
        vv = v;

        for (int k = 0; k < n; k++) {
#pragma parallel
            for (int i = 0; i < n; i++) {
#pragma ivdep
                for (int j = 0; j < n; j++) {
                    uv = (*uu) + (*vv++);

                    *xx = min(uv, *xx);
                    xx++;
                }
                vv = vv - B;
                uu = uu + B;
            }
            xx = x;
            vv = vv + B;
            uu = u + k + 1;
        }
        return;

    } else {
        int nn = (n >> 1);
        int nn2 = nn * nn;

        const int m11 = 0;
        int m12 = m11 + nn2;
        int m21 = m12 + nn2;
        int m22 = m21 + nn2;

        cilk_spawn func_B(x + m11, u + m11, v + m11, nn);
        func_B(x + m12, u + m11, v + m12, nn);
        cilk_sync;

        cilk_spawn func_D(x + m21, u + m21, v + m11, nn);
        func_D(x + m22, u + m21, v + m12, nn);
        cilk_sync;

        cilk_spawn func_B(x + m21, u + m22, v + m21, nn);
        func_B(x + m22, u + m22, v + m22, nn);

        cilk_sync;
        cilk_spawn func_D(x + m11, u + m12, v + m21, nn);
        func_D(x + m12, u + m12, v + m22, nn);
        cilk_sync;
    }
}

void func_A(TYPE *x, TYPE *u, TYPE *v, int n) {
    if (n <= B) {
        register TYPE uv;
        TYPE *xx, *uu, *vv;

        xx = x;
        uu = u;
        vv = v;

        for (int k = 0; k < n; k++) {
#pragma parallel
            for (int i = 0; i < n; i++) {
#pragma ivdep
                for (int j = 0; j < n; j++) {
                    uv = (*uu) + (*vv++);
                    *xx = min(uv, *xx);
                    xx++;
                }
                vv = vv - B;
                uu = uu + B;
            }
            xx = x;
            vv = vv + B;
            uu = u + k + 1;
        }
        return;

    } else {
        int nn = (n >> 1);
        int nn2 = nn * nn;

        const int m11 = 0;
        int m12 = m11 + nn2;
        int m21 = m12 + nn2;
        int m22 = m21 + nn2;

        func_A(x + m11, u + m11, v + m11, nn);

        cilk_spawn func_B(x + m12, u + m11, v + m12, nn);
        cilk_spawn func_C(x + m21, u + m21, v + m11, nn);
        cilk_sync;

        func_D(x + m22, u + m21, v + m12, nn);

        func_A(x + m22, u + m22, v + m22, nn);
        cilk_spawn func_B(x + m21, u + m22, v + m21, nn);
        cilk_spawn func_C(x + m12, u + m12, v + m22, nn);
        cilk_sync;
        func_D(x + m11, u + m12, v + m21, nn);
    }
}

void copy(TYPE *input, TYPE *output_1, int n)
{
    cilk_for (int i = 0; i < (n * n); ++i)
    {
        output_1[i] = input[i];
    }
    return;
}

void LOG(int * my_coords, char * msg)
{

 #ifdef _DEBUG   
    char buffer[1024];
    sprintf (buffer, "Process [%d, %d]: %s \n", my_coords[ROW], my_coords[COL], msg);
    fprintf(stderr, "%s",buffer);
#endif
}


void FW_APSP(TYPE *X, int n_local, int r, MPI_Comm comm_2d,
             MPI_Comm comm_row, MPI_Comm comm_col, int my_2d_rank, int my_coords[2])
{

    int row_root, // used for finding the source of broadcast in the row communicator(s)
        col_root, // used for finding the source of broadcast in the column communicator(s)
        coords[2]; // used for finding the source of the broadcast in the row/column communicator(s)

    TYPE *U = (TYPE *) malloc(n_local * n_local * sizeof (TYPE));
    TYPE *V = (TYPE *) malloc(n_local * n_local * sizeof (TYPE));

    int k;

    for (k = 0; k < r; ++k)
    {
        if (my_coords[ROW] == k && my_coords[COL] == k)
        {
            LOG(my_coords, "Running function A");
            // A_FW(X, n_local); // which is the shared memory function call A_FW (either serial version or the one implemented in intel cilk+)
            func_A(X, X, X, n_local);
            copy(X, U, n_local);
            copy(X, V, n_local);
        }

        if(my_coords[COL] == k)
            MPI_Bcast(V, n_local * n_local, MPI_TYPE, k, comm_col);

        if(my_coords[ROW] == k)
            MPI_Bcast(U, n_local * n_local, MPI_TYPE, k, comm_row);

        if(my_coords[ROW] == k && my_coords[COL] != k)
        {
            LOG(my_coords, "Running function B");
            // B_FW(X, U, n_local);
            func_B(X, U, V, n_local);
            copy(X, V, n_local);
        }

        if(my_coords[COL] != k)
            MPI_Bcast(V, n_local * n_local, MPI_TYPE, k, comm_col); 

        if(my_coords[ROW] != k && my_coords[COL] == k)
        {
            LOG(my_coords, "Running function C");
            // C_FW(X, V, n_local);
            func_C(X, U, V, n_local);
            copy(X, U, n_local);
        }

        if(my_coords[ROW] != k)
            MPI_Bcast(U, n_local * n_local, MPI_TYPE, k, comm_row);

        if(my_coords[1] != k  &&  my_coords[0] != k )
        {
            LOG(my_coords, "Running function D");
            // D_FW(X, U, V, n_local);
            func_D(X, U, V, n_local);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(U);
    free(V);
}



int main(int argc, char *argv[]) {
    int N_DIST = 1;

    N = 1;
    if (argc > 1)
        N_DIST = atoi(argv[1]); // Total grid size (so the grid = N_DIST X N_DIST)
    if (argc > 2)
        B = atoi(argv[2]);
    if (argc > 3) {
            if (0!= __cilkrts_set_param("nworkers",argv[3])) {
                    cout<<"Failed to set worker count\n"<<endl;
                    return 1;
            }
        }
    // cout<<"the actual worker count is "<<__cilkrts_get_nworkers()<<endl;

    int r,  // r = sqrt(num_processor)
        num_procs, // number of processes in the distributed memory setting
        dims[2], // dimensions of the cartesian process grid
        periods[2], // whether the process grid is periodic or not
        keep_dims[2]; // used for making the sub-groups (communicators for the columns and rows of the matrix)

    int my_rank, // the rank of the process in MPI_COMM_WORLD communicator in 1-d
        my_2d_rank, // the rank of the process in comm_2d communicator in 1-d
        my_coords[2]; // the coordinates of the process in 2-d process grid

    MPI_Comm comm_2d, // 2-d communicator process grid
         comm_row, // communicator which consists of the processes in the same row where the process locates
         comm_col; // communicator which consists of the processes in the same column where the process locates

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    dims[ROW] = dims[COL] = r =  static_cast<int>(std::sqrt(num_procs));
    periods[ROW] = periods[COL] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &comm_2d);

    MPI_Comm_rank(comm_2d, &my_2d_rank);
    MPI_Cart_coords(comm_2d, my_2d_rank, 2, my_coords);

    keep_dims[ROW] = 0;
    keep_dims[COL] = 1;
    MPI_Cart_sub(comm_2d, keep_dims, &comm_row);

    keep_dims[ROW] = 1;
    keep_dims[COL] = 0;
    MPI_Cart_sub(comm_2d, keep_dims, &comm_col);

    N = static_cast<int>(N_DIST / dims[ROW]);
    srand(time(0) + my_rank);

    if (B > N)
        B = N;

    NP = N;
    int NN = 2;

    while (NN < N)
        NN = NN << 1;

    if (NN == N) {
        NP = N + 1;
    }
    dist = (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

    X = (TYPE *)_mm_malloc(N * N * sizeof(TYPE), ALIGNMENT);

    // Taking average of 5 runs
    double iterGflops = 0.0;
    for (int iter = 0; iter < 5; iter++){

        cilk_for(int i = 0; i < N; i++) {
            cilk_for(int j = 0; j < N; j++) {
                if (i == j) {
                    dist[i * N + j] = 0;
                } else {
                    dist[i * N + j] = abs((j - i) % 4) + 1;
                }
            }
        }
        conv_RM_2_ZM_RM(X, 0, 0, (N), (N));
        unsigned long long tstart1 = cilk_getticks();

        long long count[1];
        int eventset=PAPI_NULL;
        int retval = 0;
        if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
            cout << "PAPI_init failed at " << my_rank << endl;
            exit(1);
        }
        retval=PAPI_create_eventset(&eventset);
        retval=PAPI_add_named_event(eventset,"PAPI_SP_OPS");
        PAPI_start(eventset);

        double t_start = MPI_Wtime();
        // funcA(X, X, X, N, N, N);
        FW_APSP(X, N, r, comm_2d, comm_row, comm_col, my_2d_rank, my_coords);
        double t_end = MPI_Wtime();

        retval=PAPI_stop(eventset,count);

        double Gflops = count[0]/(t_end - t_start);
        Gflops /= 1e9;

        // printf("\n SP: %lld %lf\n", count[0], Gflops);

        double *Gflops_all = NULL;
        if (my_rank == 0) {
            Gflops_all = (double *)malloc(sizeof(double) * num_procs);
        }
        MPI_Gather(&Gflops, 1, MPI_DOUBLE, Gflops_all, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
        
        if (my_rank == 0){
            double Gflops_sum = 0.0;
            for (int i = 0; i < num_procs; i++)
                Gflops_sum += Gflops_all[i];
            // cout << "Gflops Global: " << Gflops_sum << endl;
            iterGflops += Gflops_sum;
        }

        PAPI_shutdown();


        unsigned long long tend1 = cilk_getticks();
        conv_ZM_2_RM_RM(X, dist, 0, 0, (N), (N));
    }
    if (my_rank == 0){
        cout << num_procs << "\t" << N_DIST << "\t" << iterGflops/5.0 << endl;
    }

    // cout << "CO," << N << "," << cilk_ticks_to_seconds(tend1 - tstart1) << endl;
#ifdef pdebug
    cout << "recursive version\n";
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            cout << dist[i * N + j] << " ";
        }
        cout << "\n";
    }
#endif

    _mm_free(X);

    _mm_free(dist);

    MPI_Comm_free(&comm_col);
    MPI_Comm_free(&comm_row);
    MPI_Comm_free(&comm_2d);

    MPI_Finalize();

    return 0;
}
