//
//  path.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 02/09/2023.
//

#ifndef path_h
#define path_h

#include <vector>


unsigned long long path_hash (std::vector<unsigned long long> path, const unsigned long long NBR_STATES);

void print_path (std::vector<unsigned long long> path, bool forward=true);

#endif /* path_h */
