//
//  path.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 02/09/2023.
//

#ifndef path_h
#define path_h


unsigned long long path_hash (std::vector<unsigned long long> path, const unsigned long long NBR_STATES) {
    unsigned long long ret = 0ull;
    for (auto& state : path) {
        ret = ret*NBR_STATES + state;
    }
    return ret;
}

#endif /* path_h */
