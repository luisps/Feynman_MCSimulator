//
//  layer.hpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 18/08/2023.
//

#ifndef layer_hpp
#define layer_hpp

#include <stdio.h>
#include <random>
#include "pcg_random.hpp"

#include "circuit.h"
#include "gates.h"
#include "complex.h"
#include "myReal.h"

void layer_w (TCircuitLayer *layer, int l,
              unsigned long long current_state, unsigned long long next_state, myReal &wR, myReal &wI);

myReal layer_w_prob (TCircuitLayer *layer, int l,
              unsigned long long current_state, unsigned long long next_state, myReal &wR, myReal &wI);


myReal layer_sample (TCircuitLayer* layer, int l, unsigned long long current_state,
                    unsigned long long& next_state,
                    myReal& wR, myReal& wI,
                     pcg32& e,
                     std::uniform_real_distribution<myReal>& d, bool forwardSample=true);

#endif /* layer_hpp */
