//
//  IS_paths.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 17/08/2023.
//

#include "IS_paths.hpp"
#include <thread>
#include <vector>
#include <random>

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
    
    // thread local random number generator (seeded by a local random device)
    // see https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3551.pdf
    std::random_device rdev{};
    thread_local std::default_random_engine e{rdev()};
    std::uniform_real_distribution<float>d{0.0,1.0};  // uniform distribution in[0,1[ (float)
    
    unsigned long long s;   // sample counter
    int l;                  // layer counter
    
    float l_sumR=0.f, l_sumI=0.f;  // local summ accumulators for performance reasons
    int l_non_zero_paths = 0;      // local counter for performance reasons
    
    const int L = c->size->num_layers;   // number of layers in the circuit
    
    bool debug=false;
    
    // iteratively generate samples
    for (s=0ull; s<n_samples ; s++) {
        float path_pdf= 1.f;
        float path_wR=1.f, path_wI = 0.f;
        bool zero_power_transition = false;
        
        unsigned long long next_state=0ull, current_state = init_state;  // state before the next layer
        
        if (debug) fprintf(stderr, "Sample: %llu out of %llu\n", s, n_samples);
        
        // generate the path by stochastically sampling each layer from l=0 to l=L-2
        // the last layer (l=L-1) is handled outside the 'for' loop since it is deterministically
        // connected to 'final_state'
        for (l=0 ; l< L-1 && !zero_power_transition ; l++) {
            float wR, wI, pdf;
            
            if (debug) {
                fprintf(stderr, "\tLayer: %d out of %d\n", l, L);
                fprintf(stderr, "\tCurrent state: %llu\n", current_state);
                fprintf(stderr, "\tLayer transition...\n");
            }

            // get gate layer l
            TCircuitLayer *layer = &c->layers[l];
            
            // sample this layer for the current state,
            // returning the amplitude, pdf and next state
            pdf = layer_sample (layer, l, current_state, next_state, wR, wI, e, d);

            if (debug) {
                fprintf(stderr, "\tLayer w: %f + j %f\n", wR, wI);
                fprintf(stderr, "\tTransition p: %f\n", pdf);
                fprintf(stderr, "\tNext state: %llu\n", next_state);
            }

            if (complex_abs_square(wR, wI) <= 0.f || pdf <= 0.f) {  // I believe this should never happen
                zero_power_transition = true;
                if (debug) {
                    fprintf(stderr, "\tFinishing sample\n");
                }
            }
            else {
                path_pdf *= pdf;
                complex_multiply (path_wR, path_wI, path_wR, path_wI, wR, wI);
                if (debug) {
                    fprintf(stderr, "\tPath w up to now: %f + j %f\n", path_wR, path_wI);
                    fprintf(stderr, "\tPath Transition p: %f\n", path_pdf);
                }

            }
            
            current_state = next_state;
            if (debug) fprintf (stderr,"....\n");

        }  // layers 0 .. L-2 done
        
        // final layer (L-1)
        if (!zero_power_transition) { // path still contributes (I believer it always will)
            float wR, wI;
            
            if (debug) {
                fprintf(stderr, "\tLast Layer out of %d\n", L);
                fprintf(stderr, "\tCurrent state: %llu\n", current_state);
                fprintf(stderr, "\tFinal state: %llu\n", final_state);
                fprintf(stderr, "\tLayer deterministic transition...\n");
            }
            
            TCircuitLayer *last_layer;
            /*fprintf (stderr, "VERIFFICATION OF ALL LAYERS\n");
            for (int lll=0 ; lll<L ; lll++){
                // get gate layer L-l
                last_layer = &c->layers[lll];

                fprintf (stderr, "Layer verification: %d out of %d with %d gates", lll, L, last_layer->num_gates);
                for (int gg=0 ; gg<last_layer->num_type_gates[0] ; gg++) {
                    fprintf (stderr, "> gate %d ", (((TGate1P0 *) last_layer->gates[0])[gg]).name);
                }
                fprintf(stderr, "\n");
            }*/
            
            //fprintf (stderr, "CALLING print_circuit...\n");
            //print_circuit(c);
            //fprintf (stderr, "... print_circuit DONE...\n");

            last_layer = &c->layers[L-1];

            // evaluate this layer amplitude thansitioning from the current state to final state
            layer_w (last_layer, L-1, current_state, final_state, wR, wI);

            if (debug) {
                fprintf(stderr, "\tLast Layer w: %f + j %f\n", wR, wI);
            }

            if (complex_abs_square(wR, wI) <= 0.f) {  // This will happen frfequently
                zero_power_transition = true;
                if (debug) {
                    fprintf(stderr, "\tFinishing sample\n");
                }
            }
            else {
                complex_multiply (path_wR, path_wI, path_wR, path_wI, wR, wI);
                if (debug) {
                    fprintf(stderr, "\tPath final w : %f + j %f\n", path_wR, path_wI);
                    fprintf(stderr, "\tPath Transition p: %f\n", path_pdf);
                }
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
            if (debug) {
                fprintf(stderr, "\tPath contribution : %f + j %f\n", path_contR, path_contI);
                fprintf(stderr, "\tSamples sum: %f + j %f\n", l_sumR, l_sumI);
            }

        }
        // next sample
        if (debug) {
            fprintf(stderr, "\n");
        }

    } // end iterating over samples
    
    non_zero_paths = l_non_zero_paths;
    sumR=l_sumR;
    sumI=l_sumI;
    
    return true;
}
