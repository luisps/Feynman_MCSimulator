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

#ifdef CONVERGENCE_STATS
bool IS_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
               float &estimateR, float &estimateI, std::vector<T_Stats>& stats, const int n_threads=1);
#else
bool IS_paths (TCircuit *c, unsigned long long init_state,
                    unsigned long long final_state, const unsigned long long n_samples,
                   float &estimateR, float &estimateI, const int n_threads=1) ;
#endif


#endif /* IS_paths_hpp */
