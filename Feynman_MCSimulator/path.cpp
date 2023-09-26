//
//  path.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 25/09/2023.
//

#include <stdio.h>
#include "path.h"

void print_path (std::vector<unsigned long long> path, bool forward) {
    for (int l=(forward?0:(path.size()-1)) ; (forward ? l < path.size() : l>=0) ; l+=(forward ? 1 : -1)) {
        fprintf (stderr, "\t%llu ->", path[l]);
    }
    fprintf(stderr, "\n");
}

unsigned long long path_hash (std::vector<unsigned long long> path, const unsigned long long NBR_STATES) {
    unsigned long long ret = 0ull;
    for (auto& state : path) {
        ret = ret*NBR_STATES + state;
    }
    return ret;
}
