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
 * runtime printfs related to debugging
 *
 */

//#define DEBUG

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

#endif /* PreProcessorSettings_h */
