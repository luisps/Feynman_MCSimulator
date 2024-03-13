//
//  all_paths.hpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 31/07/2023.
//

#ifndef all_paths_hpp
#define all_paths_hpp

#include <stdio.h>
#include "circuit.h"
#include "myReal.h"

bool all_paths (TCircuit *, unsigned long long init_state,
                unsigned long long final_state,
                myReal &estimateR, myReal &estimateI, const int n_threads=1);

#endif /* all_paths_hpp */
