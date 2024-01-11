//
//  PreProcessorSettings.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 30/08/2023.
//

#ifndef PreProcessorSettings_h
#define PreProcessorSettings_h

/*
 * define the pre processor variable below to compile
 * runtime printfs related to debugging MC process
 *
 */

//#define DEBUG
//#define DEBUG_LAYER

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
    float Ar, Ai;
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
    unsigned long long n_samples;
    unsigned long long n_Paths;
    float sumR, sumI;
    long long time_us;
} T_Stats;
#endif

#endif /* PreProcessorSettings_h */
