//
//  IS_paths.hpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 17/08/2023.
//

#ifndef IS_paths_hpp
#define IS_paths_hpp

#include <stdio.h>

#include <stdio.h>
#include "circuit.h"

bool IS_paths (TCircuit *, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
                float &estimateR, float &estimateI, const int n_threads=1);


#endif /* IS_paths_hpp */
