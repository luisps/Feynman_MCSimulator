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

#include "circuit.h"
#include "gates.h"
#include "complex.h"

void layer_w (TCircuitLayer *layer, int l,
              unsigned long long current_state, unsigned long long next_state, float &wR, float &wI);

float layer_sample (TCircuitLayer* layer, int l, unsigned long long current_state,
                    unsigned long long& next_state,
                    float& wR, float& wI,
                    std::default_random_engine& e, std::uniform_real_distribution<float>& d, bool forwardSample=true);

#endif /* layer_hpp */
