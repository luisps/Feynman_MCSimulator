//
//  PreProcessorSettings.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 30/08/2023.
//

#ifndef PreProcessorSettings_h
#define PreProcessorSettings_h

#include "myReal.h"
#include "CState.h"
/*
 * define the pre processor variable below to compile
 * runtime printfs related to debugging MC process
 *
 */

//#define DEBUG
//#define DEBUG_LAYER

/*
 *  define the variable below to print a character while computing (liveness)
 */

#define __I_AM_ALIVE__

/*
 * define the pre processor variable below to compile
 * runtime printfs related to debugging the threading code
 *
 */

//#define DEBUG_THREAD


/*
 * define the pre processor variable below to compile
 * runtime code related to counting and storing non zero paths
 *
 */

#define NON_ZERO_PATHS

#ifdef NON_ZERO_PATHS
typedef struct path_desc {
    unsigned long long path_ndx;
    std::vector<unsigned long long> path;
    myReal Ar, Ai;
} T_Non_zero_path;
#endif

/*
 * define the pre processor variable below to compile
 * runtime code related to gathering statistics on convergence with nbr of samples
 *
 */

#define CONVERGENCE_STATS

#ifdef CONVERGENCE_STATS

typedef struct T_Stats_s {
    SampleCounter n_samples;
    SampleCounter n_Paths;
    myReal sumR, sumI;
    long long time_us;
#ifdef NON_ZERO_PATHS
    SampleCounter n_nonZero_paths;
#endif
} T_Stats;
#endif

#endif /* PreProcessorSettings_h */
