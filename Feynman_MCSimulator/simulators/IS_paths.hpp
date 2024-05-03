//
//  IS_paths.hpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 17/08/2023.
//

#ifndef IS_paths_hpp
#define IS_paths_hpp

#include <stdio.h>
#include <vector>

#include "circuit.h"

#include "PreProcessorSettings.h"
#include "myReal.h"
#include "CState.h"

#ifdef CONVERGENCE_STATS
bool IS_paths (TCircuit *c, CState init_state,
               CState final_state, const SampleCounter n_samples,
               myReal &estimateR, myReal &estimateI, std::vector<T_Stats>& stats, const int n_threads=1);
#else
bool IS_paths (TCircuit *c, CState init_state,
               CState final_state, const SampleCounter n_samples,
                   myReal &estimateR, myReal &estimateI, const int n_threads=1) ;
#endif


#endif /* IS_paths_hpp */
