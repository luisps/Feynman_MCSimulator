//
//  BD_paths.hpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 22/09/2023.
//

#ifndef BD_paths_hpp
#define BD_paths_hpp

#include <stdio.h>
#include <vector>

#include "circuit.h"

#include "PreProcessorSettings.h"

#ifdef CONVERGENCE_STATS
bool BD_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
               myReal &estimateR, myReal &estimateI, std::vector<T_Stats>& stats, const int n_threads=1, const bool MIS=true);
#else
bool BD_paths (TCircuit *c, unsigned long long init_state,
                    unsigned long long final_state, const unsigned long long n_samples,
                   myReal &estimateR, myReal &estimateI, const int n_threads=1, const bool MIS=true) ;
#endif




#endif /* BD_paths_hpp */
