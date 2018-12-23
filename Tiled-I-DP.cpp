#include <iostream>
#include <math.h>
#include<cilk/cilk.h>
#include <cilk/cilk_api.h>

//@Copyright: end_jsMIN Jahan Tithi, Rezaul Chowdhury, Department of Computer
// Scend_ince, Stony Brook University, Ny-11790
// Contact: jtithi@cs.stonybrook.edu,

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "cilktime.h"
#ifdef USE_PAPI
#include "papilib.h"
#endif
using namespace std;


#define MAX(x,y)    ((x) > (y)? (x) : (y))
#define MIN(x,y)    ((x) < (y)? (x) : (y))
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))

int *D, *D_loop;
int N;
int main(int argc, char *argv[]) {
    if (argc > 1)
        N = atoi(argv[1]);
    int basecase = 1024;
    if (argc > 2)
        basecase = atoi(argv[2]);
    
#ifdef USE_PAPI
    papi_init();
#endif
    /* End of CLooG code */
    if (argc > 3) {
        
        if (0 != __cilkrts_set_param("nworkers", argv[3])) {
            printf("Failed to set worker count\n");
            return 1;
        }
       
    }
    printf("#threads=%d,", __cilkrts_get_nworkers());
    D = (int *) _mm_malloc(N * N * sizeof(int), 64);
#ifdef LOOPDP
    D_loop = (int *) _mm_malloc(N * N * sizeof(int), 64);
#endif
    
    
    cilk_for(int i = 0; i < N;
             i++
             ) {
        cilk_for(int j = 0; j < N; j++) {
            if (i == j) {
                D[i*N+j] = 0;
#ifdef LOOPDP
                D_loop[i*N+j] = 0;
#endif
                
            } else {
                D[i*N+j] = abs((j - i) % 4) + 1;
#ifdef LOOPDP
                D_loop[i*N+j] =  abs((j - i) % 4) + 1;
#endif
            }
        }
    }
    
    const int baselog = log2(basecase);
    
    printf(" Block size=%d, ", basecase);
    unsigned long long tstart = cilk_getticks();
    const int end = ceil(double(N) / double(basecase));
    
    if (N >= 1) {
        int NN = 0;
        for (int k = 0; k < N; k++) {
            int *Dk = D + NN;
            cilk_for (int ii=0;ii < end;ii++) {
                const int start_i = (ii<<baselog);
                const int end_i = MIN(N, start_i + basecase);
                
                //no cilk_for
                cilk_for (int jj=0;jj<end;jj++) {
#ifdef USE_PAPI
                    int id = tid();
                    papi_for_thread(id);
                    int retval = 0;
                    if ( (retval=PAPI_start(EventSet[id])) != PAPI_OK)
                        ERROR_RETURN(retval);
#endif
                    const int start_j = (jj<<baselog);
                    const int end_j = MIN(N, start_j + basecase);
                    int nn = start_i * N;
                    
#pragma parallel
                    
                    for (int i=start_i;i<end_i;i++) {
                        int *Di = D + nn;
                        const int D4k = Di[k];
#pragma ivdep
                        for (int j=start_j;j<end_j;j++) {
                            Di[j] = MIN(Di[j], D4k + Dk[j]);
                        }
                        nn = nn + N;
                    }
#ifdef USE_PAPI
                    countMisses(id );
#endif
                }
            }
            NN = NN + N;
        }
    }
    
    unsigned long long tend = cilk_getticks();
    printf(" N=%d, Tiled runtime=%f,", N, cilk_ticks_to_seconds(tend - tstart));
    
#ifdef LOOPDP
    for(int k = 0; k < N; k++ )
    {
#pragma ivdep
        cilk_for(int i = 0; i < N; i++)
        {
            int dik = D_loop[i*N+k];
            //cilk_for
#pragma ivdep
            for(int j = 0; j < N; j++)
            {
                D_loop[i*N+j] = min(D_loop[i*N+j],  (dik + D_loop[k*N+j]));
                
            }
        }
    }
    for(int i =0; i<N;i++)
    {
        for(int j = 0; j <N; j++){
            assert(D_loop[i*N+j]==D[i*N+j]);
        }
    }
#endif
    
#ifdef USE_PAPI
    countTotalMiss(p);
    PAPI_shutdown();
    delete threadcounter;
    for (int i = 0; i < p; i++) delete l2miss[i];
    delete l2miss;
    delete errstring;
    delete EventSet;
    delete eventCode;
#endif
#ifdef ADAPTIVITY
    cout << "," << tstart << " ," << tend << endl;
#endif
    _mm_free(D);
    cout << endl;
    return 0;
}
