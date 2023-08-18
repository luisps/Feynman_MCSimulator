//
//  IS_paths.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 17/08/2023.
//

#include "IS_paths.hpp"
#include <thread>
#include <vector>

#include "complex.h"
#include "layer.hpp"

static bool IS_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                float &sumR, float &sumI, int& non_zero_paths);


bool IS_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state, const unsigned long long n_samples,
               float &estimateR, float &estimateI, const int n_threads) {
    
    bool ret=true;
    int non_zero_paths = 0;
    float sumR=0.f, sumI=0.f;

    
    if (n_threads<=1) {
        ret = IS_paths_NOVEC_T (c, init_state, final_state, n_samples, sumR, sumI, non_zero_paths);
        estimateR = sumR / ((float)n_samples);
        estimateI = sumI / ((float)n_samples);
    }
    else {
        const unsigned long long ln_samples = (const unsigned long long)(n_samples / (unsigned long long)  n_threads);
        if (ln_samples <1) {
            fprintf (stderr, "Error: some threads would be assigned 0 samples!\nABORT!");
            ret = false;
        }
        else
        {
            std::vector<std::thread> threads;
            float * l_sumR = new float[n_threads];
            float * l_sumI = new float [n_threads];
            int * l_NzeroP = new int [n_threads];

            for (int t = 0; t < n_threads; t++) {
                l_sumR[t] = l_sumI[t] = 0.f;
                l_NzeroP[t] = 0;
                threads.push_back(std::thread(IS_paths_NOVEC_T, c, init_state, final_state, ln_samples, std::ref(l_sumR[t]), std::ref(l_sumI[t]), std::ref(l_NzeroP[t])));
            }
            
            for (auto &th : threads) {
                th.join();
            }
            
            for (int t=0 ; t< n_threads ; t++) {
                non_zero_paths += l_NzeroP[t];
                sumR += l_sumR[t];
                sumI += l_sumI[t];
            }
            
            estimateR = sumR / n_samples;
            estimateI = sumI / n_samples;
            
            delete[] l_NzeroP;
            delete[] l_sumI;
            delete[] l_sumR;
        }
    }

    fprintf (stderr, "Non zero paths: %d\n", non_zero_paths);

    return ret;
}


static bool IS_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state, const unsigned long long n_samples,
                              float &sumR, float &sumI, int& non_zero_paths) {
    
    unsigned long long s;   // sample counter
    int l;                  // layer counter
    
    float l_sumR=0.f, l_sumI=0.f;  // local summ accumulators for performance reasons
    int l_non_zero_paths = 0;      // local counter for performance reasons
    
    const int L = c->size->num_layers;   // number of layers in the circuit
    
    // iteratively generate samples
    for (s=0ull; s<n_samples ; s++) {
        float path_pdf= 1.f;
        float path_wR=1.f, path_wI = 0.f;
        bool zero_power_transition = false;
        
        unsigned long long next_state=0ull, current_state = init_state;  // state before the next layer
        
        // generate the path by stochastically sampling each layer from l=0 to l=L-2
        // the last layer (l=L-1) is handled outside the 'for' loop since it is deterministically
        // connected to 'final_state'
        for (l=0 ; l< L-2 && !zero_power_transition ; l++) {
            float wR, wI, pdf;
            
            // get gate layer l
            TCircuitLayer *layer = &c->layers[l];
            
            // sample this layer for the current state,
            // returning the amplitude, pdf and next state
            pdf = layer_sample (layer, l, current_state, next_state, wR, wI);
            
            if (complex_abs_square(wR, wI) <= 0.f || pdf <= 0.f) {  // I believe this should never happen
                zero_power_transition = true;
            }
            else {
                path_pdf *= pdf;
                complex_multiply (path_wR, path_wI, path_wR, path_wI, wR, wI);
            }
            
            current_state = next_state;

        }  // layers 0 .. L-2 done
        
        // final layer (L-1)
        if (!zero_power_transition) { // path still contributes (I believer it always will)
            float wR, wI;
            
            // get gate layer L-l
            TCircuitLayer *layer = &c->layers[L-l];
            
            // evaluate this layer amplitude thansitioning from the current state to final state
            layer_w (layer, L-1, current_state, final_state, wR, wI);

            if (complex_abs_square(wR, wI) <= 0.f) {  // This will happen frfequently
                zero_power_transition = true;
            }
            else {
                complex_multiply (path_wR, path_wI, path_wR, path_wI, wR, wI);
            }

        }
        
        if (!zero_power_transition) {  // OK, count non_zero paths and accumulate
            l_non_zero_paths++;
            float pdf_reciprocal = 1.f / path_pdf;
            float path_contR, path_contI;
            complex_multiply (path_contR, path_contI, path_wR, path_wI, pdf_reciprocal, 0.f);
            
            // accumulate on local sums
            l_sumR += path_contR;
            l_sumI += path_contI;
            
        }
    } // end iterating over samples
    
    non_zero_paths = l_non_zero_paths;
    sumR=l_sumR;
    sumI=l_sumI;
    
    return true;
}
