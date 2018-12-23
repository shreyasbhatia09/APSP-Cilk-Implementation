#include <iostream>
#include <papi.h>

using namespace std;
int p = 0;
#define NUM_EVENTS 3

#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); }
bool *threadcounter;
long long int **l2miss;
double *gflops;
char **errstring;
int *EventSet, *eventCode;



inline unsigned long tid() {
    return __cilkrts_get_worker_number();
}


void inline papi_for_thread(int id) {
    int retval = 0;
    if (threadcounter[id] == false) {
        if (PAPI_thread_init(tid) != PAPI_OK  || PAPI_thread_id() == (unsigned long int) -1) {
            fprintf(stderr, "Error: %s\n", errstring[id]);
            exit(1);
        }
        /* Creating event set   */
        EventSet[id] = PAPI_NULL;

        retval=PAPI_create_eventset(&EventSet[id]);
        if(retval != PAPI_OK)
        {
            cout<<"cannot create papi event"<<endl;
            ERROR_RETURN(retval);
        }
        retval = PAPI_add_event(EventSet[id], PAPI_SP_OPS);
//        retval = PAPI_add_named_event(EventSet[id],"PAPI_SP_OPS");
        if ( ( retval ) != PAPI_OK )
        {
            cout<<"Add event failed"<<endl;
            ERROR_RETURN(retval);

        }

        /* Add the array of events PAPI_TOT_INS and PAPI_TOT_CYC to the eventset*/
//        if ((retval = PAPI_add_event(EventSet[id], PAPI_L1_TCM)) != PAPI_OK)
//            ERROR_RETURN(retval);

        //if ((retval = PAPI_add_event(EventSet[id], PAPI_L1_ICM)) != PAPI_OK)
        //    ERROR_RETURN(retval);
//
//        if ((retval = PAPI_add_event(EventSet[id], PAPI_L2_TCM)) != PAPI_OK)
//        ERROR_RETURN(retval);

        //if ((retval = PAPI_add_event(EventSet[id], PAPI_L2_ICM)) != PAPI_OK)
        //    ERROR_RETURN(retval);

//        if ((retval = PAPI_add_event(EventSet[id], PAPI_L3_TCM)) != PAPI_OK)
//        ERROR_RETURN(retval);

        /* Start counting */


        threadcounter[id] = true;

    }
}

void inline countMisses(int id) {
//    long long int comiss[NUM_EVENTS];
//    int retval = 0;
//    if ((retval = PAPI_stop(EventSet[id], comiss)) != PAPI_OK) ERROR_RETURN(retval);
//
//    for (int i = 0; i < NUM_EVENTS; i++) {
//        l2miss[id][i] += comiss[i];
//    }
}

void inline countGflops(int id)
{
    double iterGflops = 0.0;
    long long count[1];
    int retval = 0;

    if ((retval = PAPI_stop(EventSet[id], count)) != PAPI_OK) ERROR_RETURN(retval);
    double Gflops = count[0];
    Gflops /= 1e9;
    gflops[id] += Gflops;
    //printf("\n SP: %lld %lf\n", count[0], Gflops);
    return;
}

void papi_reset() {
    p = __cilkrts_get_nworkers();
    for (int i = 0; i < p; i++) {
        PAPI_reset(EventSet[i]);
        l2miss[i][0:NUM_EVENTS]=0;
    }
}

void inline papi_init() {

    p = __cilkrts_get_nworkers();
    EventSet = new int[p];
    EventSet[0:p]=PAPI_NULL;

    for (int i = 0; i < p; i++)
        PAPI_reset(EventSet[i]);
    l2miss = new long long int *[p];
    threadcounter = new bool[p];
    threadcounter[0:p]=false;

    errstring = new char *[p];

    gflops = new double[p];

    eventCode = new int[p];
    eventCode[0:p]=0;
    for (int i = 0; i < p; i++) {
        l2miss[i] = new long long int[NUM_EVENTS];

        l2miss[i][0:NUM_EVENTS]=0;
        errstring[i] = new char[PAPI_MAX_STR_LEN];
        gflops[i] = 0.0;
    }

    int id = __cilkrts_get_worker_number();
    int retval = 0;
    if ((retval = PAPI_library_init(PAPI_VER_CURRENT)) != PAPI_VER_CURRENT) {
        fprintf(stderr, "Error: %s\n", errstring[id]);
        exit(1);
    }

}

double countTotalGflops(int p) {

    double total = 0.0;
    for (int i = 0; i < p; i++)
    {
        total += gflops[i];
    }
    return total;
}
void countTotalMiss(int p) {
    long long int l1cmiss_total[NUM_EVENTS];
    l1cmiss_total[0:NUM_EVENTS]=0;
    int retval = 0;
    for (int j = 0; j < NUM_EVENTS; j++) {
        for (int i = 0; i < p; i++) {


            l1cmiss_total[j] += l2miss[i][j];

            if (threadcounter[i] == true) {
                if ((retval = PAPI_remove_event(EventSet[i], PAPI_L1_TCM))
                    != PAPI_OK) ERROR_RETURN(retval);

                /* if ((retval = PAPI_remove_event(EventSet[i], PAPI_L1_ICM))
                     != PAPI_OK)
                     ERROR_RETURN(retval);*/

                if ((retval = PAPI_remove_event(EventSet[i], PAPI_L2_TCM))
                    != PAPI_OK) ERROR_RETURN(retval);

                /*if ((retval = PAPI_remove_event(EventSet[i], PAPI_L2_ICM))
                    != PAPI_OK)
                    ERROR_RETURN(retval);*/

                if ((retval = PAPI_remove_event(EventSet[i], PAPI_L3_TCM))
                    != PAPI_OK) ERROR_RETURN(retval);

                /* Free all memory and data structures, EventSet must be empty. */
                if ((retval = PAPI_destroy_eventset(&EventSet[i])) != PAPI_OK) ERROR_RETURN(retval);

                //PAPI_L3_TCM
                /* free the resources used by PAPI */

                threadcounter[i] = false;
            }
        }
        cout << "L" << j + 1 << " TCM = " << l1cmiss_total[j] << ", ";
    }

}

void papi_end() {

    delete[]threadcounter;
    for (int i = 0; i < p; i++) {
        delete[]l2miss[i];
        delete[]errstring[i];
    }
    delete[]EventSet;
    delete[]eventCode;
    PAPI_shutdown();
}
///////////////////////////////////////// utility functions for papi ends here//////////////////////////////
