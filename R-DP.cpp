//@Copyright: JesMIN Jahan Tithi, Rezaul Chowdhury, Department of Computer Science, Stony Brook University, Ny-11790
//Contact: jtithi@cs.stonybrook.edu,

//Aakrati Talati : Modified to accept any N (Row Major layout) May, 2013

//Optimized: JesMIN Jahan Tithi, Date: Dec 20, 2013

//compile with :icpc -O3 -AVX  -xhost -o  fwr FW_rec.cpp  -ansi-alias -ip -opt-subscript-in-range -opt-prefetch=4 -fomit-frame-pointer -funroll-all-loops -vec-report  -parallel -restrict


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include<cilk/cilk.h>
#include <cilk/cilk_api.h>

#include <iostream>

#include "cilktime.h"
using namespace std;

#ifdef USE_PAPI
#include "papilib.h"
#endif

#ifdef ENERGY
#include "energylib.h"
#endif

#ifndef TYPE
#define TYPE float
#endif

#ifndef ALIGNMENT
#define ALIGNMENT 64
#endif

int N, B, NP, NPDQ;
TYPE *dist, *X, *D;

#define MIN(a, b) (a<b?a:b)
#define MAX(a, b) (a>b?a:b)


void funcD( TYPE *x, TYPE *u, TYPE *v, int n, int xi, int xj, int uj) {
    if (xi >= N || xj >= N || uj >= N)
        return;

    if (n <= B) {
#ifdef USE_PAPI
        int id = tid();
        int retval = 0;
        papi_for_thread(id);
        if ( (retval=PAPI_start(EventSet[id])) != PAPI_OK)
            ERROR_RETURN(retval);
#endif
        int endi = (xi + n >= N) ? (N - xi) : n;
        int endj = (xj + n >= N) ? (N - xj) : n;
        int endu = (uj + n >= N) ? (N - uj) : n;

        // Copy Optimization.
        __declspec(align(64)) TYPE V[endu * endj];
//#pragma parallel
//        for (int i = 0; i < endu; i++) {
//#pragma ivdep
//            for (int j = 0; j < endj; j++) {
//                V[j * endu + i] = v[i * endj + j];
//            }
//        }
        for (int k = 0; k < endu; k++) {
#pragma parallel
            for (int i = 0; i < endi; i++) {
#pragma ivdep
                for (int j = 0; j < endj; j++) {
                    x[i*NPDQ+j] = MIN((u[i*NPDQ+k]) + (v[k*NP+j]), x[i*NPDQ+j]);
                }
                x[i*NPDQ + j] = x_ij;
            }
        }
#ifdef USE_PAPI
        countGflops(id);
//        countMisses(id );
#endif
        return;
    } else {
        int nn = (n >> 1);

        int nn2 = nn * nn;

        const int m11 = 0;
        int m12 = m11 + nn;
        int m21 = nn * NPDQ;
        int m22 = m21 + nn;

        cilk_spawn funcD(x + m11, u + m11, v + m11, nn, xi, xj, uj);
        cilk_spawn funcD(x + m12, u + m11, v + m12, nn, xi, xj + nn, uj);
        cilk_spawn funcD(x + m21, u + m21, v + m11, nn, xi + nn, xj, uj);
        funcD(x + m22, u + m21, v + m12, nn, xi + nn, xj + nn, uj);

        cilk_sync;

        cilk_spawn funcD(x + m11, u + m12, v + m21, nn, xi, xj, uj + nn);
        cilk_spawn funcD(x + m12, u + m12, v + m22, nn, xi, xj + nn, uj + nn);
        cilk_spawn funcD(x + m21, u + m22, v + m21, nn, xi + nn, xj, uj + nn);
        funcD(x + m22, u + m22, v + m22, nn, xi + nn, xj + nn, uj + nn);

        cilk_sync;
    }
}

void funcC( TYPE *x, TYPE *u, TYPE *v, int n, int xi, int xj, int uj) {
    if (xi >= N || xj >= N || uj >= N)
        return;

    if (n <= B) {
#ifdef USE_PAPI
        int id = tid();
        int retval = 0;
        papi_for_thread(id);
        if ( (retval=PAPI_start(EventSet[id])) != PAPI_OK)
            ERROR_RETURN(retval);
#endif
        int endi = (xi + n >= N) ? (N - xi) : n;
        int endj = (xj + n >= N) ? (N - xj) : n;
        int endu = (uj + n >= N) ? (N - uj) : n;

        for (int k = 0; k < endu; k++) {
#pragma parallel
            for (int i = 0; i < endi; i++) {
#pragma ivdep
                for (int j = 0; j < endj; j++) {
                    x[i*NPDQ+j] = MIN((u[i*NPDQ+k]) + (v[k*NP+j]), x[i*NPDQ+j]);
                }
            }
        }
#ifdef USE_PAPI
        countGflops(id);
        //countMisses(id );
#endif
        return;
    } else {
        int nn = (n >> 1);
        int nn2 = nn * nn;

        const int m11 = 0;
        int m12 = m11 + nn;
        int m21 = nn * NPDQ;
        int m22 = m21 + nn;

        cilk_spawn funcC(x + m11, u + m11, v + m11, nn, xi, xj, uj);
        funcC(x + m21, u + m21, v + m11, nn, xi + nn, xj, uj);
        cilk_sync;
        cilk_spawn funcD(x + m12, u + m11, v + m12, nn, xi, xj + nn, uj);
        funcD(x + m22, u + m21, v + m12, nn, xi + nn, xj + nn, uj);
        cilk_sync;

        cilk_spawn funcC(x + m12, u + m12, v + m22, nn, xi, xj + nn, uj + nn);
        funcC(x + m22, u + m22, v + m22, nn, xi + nn, xj + nn, uj + nn);
        cilk_sync;

        cilk_spawn funcD(x + m11, u + m12, v + m21, nn, xi, xj, uj + nn);
        funcD(x + m21, u + m22, v + m21, nn, xi + nn, xj, uj + nn);
        cilk_sync;
    }
}

void funcB( TYPE *x, TYPE *u, TYPE *v, int n, int xi, int xj, int uj) {
    if (xi >= N || xj >= N || uj >= N)
        return;

    if (n <= B) {
#ifdef USE_PAPI
        int id = tid();
        int retval = 0;
        papi_for_thread(id);
        retval = PAPI_start(EventSet[id]);
        if ( retval != PAPI_OK)
        {
            ERROR_RETURN(retval);
        }

#endif

        int endi = (xi + n >= N) ? (N - xi) : n;
        int endj = (xj + n >= N) ? (N - xj) : n;
        int endu = (uj + n >= N) ? (N - uj) : n;

        for (int k = 0; k < endu; k++) {
#pragma parallel
            for (int i = 0; i < endi; i++) {
#pragma ivdep
                for (int j = 0; j < endj; j++) {
                    x[i*NPDQ+j] = MIN(x[i*NPDQ+j], (u[i*NPDQ+k]) + (v[k*NP+j]));
                }
            }
        }
#ifdef USE_PAPI
        countGflops(id );
//        countMisses(id );
#endif
        return;
    } else {
        int nn = (n >> 1);
        int nn2 = nn * nn;

        const int m11 = 0;
        int m12 = m11 + nn;
        int m21 = nn * NPDQ;
        int m22 = m21 + nn;

        cilk_spawn funcB(x + m11, u + m11, v + m11, nn, xi, xj, uj);
        funcB(x + m12, u + m11, v + m12, nn, xi, xj + nn, uj);
        cilk_sync;

        cilk_spawn funcD(x + m21, u + m21, v + m11, nn, xi + nn, xj, uj);
        funcD(x + m22, u + m21, v + m12, nn, xi + nn, xj + nn, uj);
        cilk_sync;

        cilk_spawn funcB(x + m21, u + m22, v + m21, nn, xi + nn, xj, uj + nn);
        funcB(x + m22, u + m22, v + m22, nn, xi + nn, xj + nn, uj + nn);
        cilk_sync;
        cilk_spawn funcD(x + m11, u + m12, v + m21, nn, xi, xj, uj + nn);
        funcD(x + m12, u + m12, v + m22, nn, xi, xj + nn, uj + nn);
        cilk_sync;
    }
}

void funcA( TYPE *x, TYPE *u, TYPE *v, int n, int xi, int xj, int uj) {
    if (xi >= N || xj >= N || uj >= N)
        return;

    if (n <= B) {
#ifdef USE_PAPI
        int id = tid();
        int retval = 0;
        papi_for_thread(id);
        if ( ( retval=PAPI_start(EventSet[id] ) ) != PAPI_OK)
            ERROR_RETURN(retval);
#endif

        int endi = (xi + n >= N) ? (N - xi) : n;
        int endj = (xj + n >= N) ? (N - xj) : n;
        int endu = (uj + n >= N) ? (N - uj) : n;

        for (int k = 0; k < endu; k++) {
#pragma parallel
            for (int i = 0; i < endi; i++) {
#pragma ivdep
                for (int j = 0; j < endj; j++) {
                    x[i*NPDQ+j] = MIN(x[i*NPDQ+j], (u[i*NPDQ+k]) + (v[k*NP+j]));

                }
            }
        }
#ifdef USE_PAPI
        countGflops(id);
//        countMisses(id );
#endif
        return;

    } else {
        int nn = (n >> 1);
        int nn2 = nn * nn;

        const int m11 = 0;
        int m12 = m11 + nn;
        int m21 = nn * NPDQ;
        int m22 = m21 + nn;

        funcA(x + m11, u + m11, v + m11, nn, xi, xj, uj);

        cilk_spawn funcB(x + m12, u + m11, v + m12, nn, xi, xj + nn, uj);
        funcC(x + m21, u + m21, v + m11, nn, xi + nn, xj, uj);
        cilk_sync;

        funcD(x + m22, u + m21, v + m12, nn, xi + nn, xj + nn, uj);

        funcA(x + m22, u + m22, v + m22, nn, xi + nn, xj + nn, uj + nn);
        cilk_spawn funcB(x + m21, u + m22, v + m21, nn, xi + nn, xj, uj + nn);
        funcC(x + m12, u + m12, v + m22, nn, xi, xj + nn, uj + nn);
        cilk_sync;
        funcD(x + m11, u + m12, v + m21, nn, xi, xj, uj + nn);
    }
}

#ifdef LOOPDP

void FloydWarshall(int n) {
    for (int k = 0; k < n; k++) {
        int ink = k * NP;
        // #pragma cilk grainsize = 256
        cilk_for(int i = 0; i < n; i++) {
#ifdef USE_PAPI
            int id = tid();
            int retval = 0;
            papi_for_thread(id);
            if ((retval = PAPI_start(EventSet[id])) != PAPI_OK) ERROR_RETURN(retval);
#endif
            int in = i * NP;
            int *d = D + in;
            int v = D[in + k];
            int *dk = D + ink;

#pragma ivdep
            for (int j = 0; j < n; j++) {
                {
                    *d = min(*d, v + *dk);
                    dk++;
                    d++;
                }
            }
#ifdef USE_PAPI
            countGflops(id);
            //countMisses(id);
#endif
        }
    }
}
#endif
int main(int argc, char *argv[]) {
    N = 0;
    B = 0;
    if (argc > 1)
        N = atoi(argv[1]);
    if (argc > 2)
        B = atoi(argv[2]);

    if (argc > 3) {

        if (0 != __cilkrts_set_param("nworkers", argv[3])) {
            printf("Failed to set worker count\n");
            return 1;
        }

    };
    printf("#threads=%d,", __cilkrts_get_nworkers());

    NP = NPDQ = N;
    int NN = 2;

    while (NN < N)
        NN = NN << 1;


#ifdef USE_PAPI
    papi_init();
#endif

    if (B > N)
        B = N;

    printf(" Basecase=%d, ", B);
#ifdef CO
    X = ( TYPE * ) _mm_malloc( NPDQ * NPDQ * sizeof( TYPE ), ALIGNMENT );
#endif

#ifdef LOOPDP
    D = (TYPE *) _mm_malloc( NP * NP * sizeof( TYPE ), ALIGNMENT );
#endif

    cilk_for(int i = 0; i < N;
    i++
    )
    {
        cilk_for(int j = 0; j< N; j++ )
        {

            if( i == j )
            {
#ifdef CO
                X[i*NPDQ+j] = 0;
#endif
#ifdef LOOPDP
                D[i*NP+j]=0;
#endif
            }
            else
            {
#ifdef CO
                X[i*NPDQ+j]=abs( (j-i)%4 + 1 ) % 1;
#endif
#ifdef LOOPDP
                D[i*NP+j]=abs( (j-i)%4 + 1 ) % 1;;
#endif
            }
        }
    }
    unsigned long long tstart, tend;
#ifdef ENERGY
    energy_lib_init();
#endif

#ifdef CO
    #ifdef ENERGY
    start_energy_measurement();
#endif
    tstart = cilk_getticks();

    funcA(X, X, X, NN,0,0,0);

    tend = cilk_getticks();
    cout<<" CO: N="<<N<<", runtime="<<cilk_ticks_to_seconds(tend-tstart)<<",";
#ifdef ENERGY
    stop_energy_measurement();
#endif
#ifdef pdebug
    cout<<"loop version\n";
    for(int i=0;i<N;i++) {
        for(int j=0;j<N;j++) {
            cout<<D[i*NP+j]<<" ";
        }
        cout<<"\n";
    }
#endif
#endif
#ifdef LOOPDP

    tstart = cilk_getticks();
    FloydWarshall(N);
    tend = cilk_getticks();

    cout<<" LOOPDP: N="<<N<<", runtime="<<cilk_ticks_to_seconds(tend-tstart);

#ifdef pdebug
    cout<<"recursive version\n";
    for(int i=0;i<N;i++) {
        for(int j=0;j<N;j++) {
            cout<<X[i*NPDQ+j]<<" ";
        }
        cout<<"\n";
    }
#endif

#ifdef CO
    cout<<",";
    for(int i = 0; i <N; i++ )
    {
        for(int j = 0; j<N;j++ )	//
        {
            if(D[i*NP+j]!=X[i*NPDQ+j])
            {
                cout<<"Wrong at"<<i<<j<<endl;
                break;
            }

        }
    }
#endif
    _mm_free(D);

#endif

#ifdef CO
    _mm_free(X);
#endif

#ifdef USE_PAPI
    double totalflops = countTotalGflops(p);
    cout<<"Total Gflops: "<<totalflops<<endl;
    cout<<"GLops per sec: "<<totalflops/cilk_ticks_to_seconds(tend-tstart)<<endl;
    cout<<"Number of Cilk Ticks: "<<tend-tstart<<endl;
    //countTotalMiss(p);
    PAPI_shutdown();
    delete threadcounter;
    for (int i = 0; i < p; i++) delete l2miss[i];
    delete l2miss;
    delete errstring;
    delete EventSet;
    delete eventCode;
#endif
    cout<<endl;

#ifdef ENERGY
    clear_energy_measurement();
#endif
    return 0;
}
