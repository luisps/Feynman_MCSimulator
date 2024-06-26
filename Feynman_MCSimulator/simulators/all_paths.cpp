//
//  all_paths.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 31/07/2023.
//

#include <stdlib.h>
#include <thread>

#include <vector>

using namespace std;

#include "all_paths.hpp"
#include "gates.h"
#include "complex.h"

#include "layer.hpp"
#include "PreProcessorSettings.h"
#include "path.h"

#ifdef NON_ZERO_PATHS
static bool all_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state,
                const unsigned long long l0_init, const unsigned long long l0_finish,
                myReal &estimateR, myReal &estimateI, SampleCounter& non_zero_paths);
#else
static bool all_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state,
                const unsigned long long l0_init, const unsigned long long l0_finish,
                myReal &estimateR, myReal &estimateI);
#endif

bool all_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state,
                myReal &estimateR, myReal &estimateI, const int n_threads) {
    
    bool ret=true;
#ifdef NON_ZERO_PATHS
    SampleCounter non_zero_paths = 0;
#endif
    // compute NBR_STATES : the number of states to iterate in each intermediate state layer
    const unsigned long long NBR_STATES = 1 << c->size->num_qubits;
    

    // all_paths currently requires at least 4 layers
    // terminate otherwise
    if (c->size->num_layers<4) {
        fprintf (stderr, "all_paths currently requires at least 4 layers!\n");
        fprintf (stderr, "The specified circuit has %d < 4 layers! Terminating!\n\n", c->size->num_layers);
        estimateR = estimateI =0.f;
        return false;
    }
    
    
    if (n_threads<=1) {
#ifdef NON_ZERO_PATHS
        ret = all_paths_NOVEC_T (c, init_state, final_state, 0, NBR_STATES, estimateR, estimateI, non_zero_paths);
#else
        ret = all_paths_NOVEC_T (c, init_state, final_state, 0, NBR_STATES, estimateR, estimateI);
#endif
    }
    else {
        const SampleCounter N_ROWS = NBR_STATES / n_threads;
        if (N_ROWS <1) {
            fprintf (stderr, "Error: some threads would be assigned 0 states at level 0!\nABORT!");
            ret = false;
        }
        else
        {
            std::vector<std::thread> threads;
            myReal * l_estimateR = new myReal[n_threads];
            myReal * l_estimateI = new myReal [n_threads];
#ifdef NON_ZERO_PATHS
            SampleCounter * l_NzeroP = new SampleCounter [n_threads];
#endif
            for (int t = 0; t < n_threads; t++) {
                l_estimateR[t] = l_estimateI[t] = 0.f;
#ifdef NON_ZERO_PATHS
                l_NzeroP[t] = 0;
                threads.push_back(std::thread(all_paths_NOVEC_T, c, init_state, final_state, N_ROWS*t, N_ROWS*(t+1), std::ref(l_estimateR[t]), std::ref(l_estimateI[t]), std::ref(l_NzeroP[t])));
#else
                threads.push_back(std::thread(all_paths_NOVEC_T, c, init_state, final_state, N_ROWS*t, N_ROWS*(t+1), std::ref(l_estimateR[t]), std::ref(l_estimateI[t])));
#endif
            }
            
            for (auto &th : threads) {
                th.join();
            }
            
            estimateR = estimateI = 0.f;
#ifdef NON_ZERO_PATHS
            non_zero_paths = 0;
#endif
            for (int t=0 ; t< n_threads ; t++) {
#ifdef NON_ZERO_PATHS
                non_zero_paths += l_NzeroP[t];
#endif
                estimateR += l_estimateR[t];
                estimateI += l_estimateI[t];
            }
            
#ifdef NON_ZERO_PATHS
            delete[] l_NzeroP;
#endif
            delete[] l_estimateI;
            delete[] l_estimateR;
        }
    }

#ifdef NON_ZERO_PATHS
    fprintf (stderr, "Non zero paths: %llu\n", non_zero_paths);
#endif
    return ret;
}

// multithreaded ready all_paths version
#ifdef NON_ZERO_PATHS
static bool all_paths_NOVEC_T (TCircuit *c,
                unsigned long long init_state, unsigned long long final_state,
                const unsigned long long l0_init, const unsigned long long l0_finish,
                myReal &estimateR, myReal &estimateI, SampleCounter& non_zero_paths) {
#else
    static bool all_paths_NOVEC_T (TCircuit *c,
                                   unsigned long long init_state, unsigned long long final_state,
                    const unsigned long long l0_init, const unsigned long long l0_finish,
                    myReal &estimateR, myReal &estimateI) {
#endif
#ifdef NON_ZERO_PATHS
    non_zero_paths = 0;
#endif
    myReal sumR=0., sumI=0.;
    const int L = c->size->num_layers;  // number of gate layers
    /* the number of intermediate states layers is L-1.
       This is the number of 'for' loops we have to realize
       We will do L-3 for loops in a single aggregated while loop:
       . the outer most loop (l=0) is implemented as a for loop,
         thus allowing FOR LOOP parallelism
         the counter associated with this outer 'for' loop is ndxs[0]
       . the inner most 'for' loop (corresponding to L-1 (just before transitioning
         to the final state) is still implememnted with a 'for' loop
         the counter associated with this inner most 'for' is 'next_state'
     */
    
    /* Note that the loop counter for each intermediate state layer encodes
     the state itself. */
    
    unsigned long long *ndxs = new unsigned long long [L-2];

    // store the intermediate amplitudes associated with each transition between layers
    myReal *w_intermediateR = new myReal [L-2];
    myReal *w_intermediateI = new myReal [L-2];

    // compute NBR_STATES : the number of states to iterate in each intermediate state layer
    const unsigned long long NBR_STATES = 1 << c->size->num_qubits;
    
    
    // outer most for loop
    for (ndxs[0]=l0_init; ndxs[0]<l0_finish ; ndxs[0]++) {
        
        // initialize all our L-3 loop counters (excludes ndxs[0])
        for (int l=1 ; l< (L-2) ; l++) ndxs[l] = 0;
        
        // We will enter the multiple variable number of for loops below
        // We must first make all initializations corresponding to all ndxs[>0]==0
        // Note that this corresponds to computing the amplitude of the path
        // transitioning from the initial state, through all intermediate states |ndxs[0]> |0>
        // all way to the final state.
        // However, since the final transition to the final state is supported by an
        // explicit 'for' loop we don't pre compute the very last one
        
        unsigned long long current_state = init_state;
        
        // Handle layer 0 explicitly to allow for ealy termination ('break')
        {
            myReal wR, wI;
            
            // get gate layer l
            TCircuitLayer *layer = &c->layers[0];
            
            // get the amplitude of this layer for the given current and next states
            layer_w (layer, 0, current_state, ndxs[0], wR, wI);
            
            if (complex_abs_square(wR, wI) <= 0.) {
                ///// VERIFY THIS
                continue;
            }

            
#if defined(DEBUG) && !defined(__FLOAT_MP__)
            fprintf (stderr, "Layer l=0 -> %f + i %f\n", wR, wI);
#endif
            
            w_intermediateR[0] = wR;
            w_intermediateI[0] = wI;
            
            // update current state
            current_state = ndxs [0];
        }  // end layer 0
        
        for (int l=1 ; l < (L-2) ; l++) {
            myReal wR, wI;
            
            // get gate layer l
            TCircuitLayer *layer = &c->layers[l];
            
            // get the amplitude of this layer for the given current and next states
            layer_w (layer, l, current_state, ndxs[l], wR, wI);
            
#if defined(DEBUG) && !defined(__FLOAT_MP__)
            fprintf (stderr, "Layer l=%d -> %f + i %f\n", l, wR, wI);
#endif
            
            w_intermediateR[l] = wR;
            w_intermediateI[l] = wI;
            
            // update current state
            current_state = ndxs [l];
        }  // end iterating over layers 1 .. L-3 (that is L-1 layers)
        
        //the code corresponding to the first iteration of outer nested L-2 for loops
         //has been executed
         // Now enter the loops , starting with the inner most for loop, i.e.,
         //from current_state ndxs[L-3] ->  all s in last intermediate state and then to final_state  \=>
         
        
        bool stop = false;
        while (!stop) {
            myReal acc_wR, acc_wI;
            acc_wR = w_intermediateR[0]; acc_wI=w_intermediateI[0];
            
#if defined(DEBUG)
            // DEBUG
            for (int l=0 ; l<L-2 ; l++) fprintf (stderr, "qb[%d]=%llu ", l, ndxs[l]);
            fprintf(stderr, "\n");
            fflush(stderr);
#endif
            
            // check whether the amplitude from previous layers is 0
            for (int l=1 ; l<(L-2) ; l++) {
                complex_multiply (acc_wR, acc_wI, acc_wR, acc_wI, w_intermediateR[l], w_intermediateI[l]);
            }
            
            // check whether accumulated amplitude is zero
            if (complex_abs_square(acc_wR, acc_wI) > 0.) {
                // propagation of current state through layer L-2
                // and immediately afterwards to final state
                // these are 2 final hops of a path
                // next state is iterated here over all possible states
                
                // layer L-2
                const int l = L-2;
                // get gate layer L-2 and lastLayer
                TCircuitLayer *layer_2L = &c->layers[l];
                TCircuitLayer *layer_last = &c->layers[l+1];
                
                current_state = (L >= 4 ? ndxs[L-3] : ndxs[0]);
                for (unsigned long long next_state=0; next_state<NBR_STATES ; next_state++) {
                    myReal path_wR, path_wI;
                    myReal wRl, wIl, wRll, wIll;
                    
                    layer_w(layer_2L, l, current_state, next_state, wRl, wIl);
#if defined(DEBUG) && !defined(__FLOAT_MP__)
                    fprintf (stderr, "Layer l=%d -> %f + i %f\n", l, wRl, wIl);
#endif
                    // check whether this layer amplitude is zero
                    if (complex_abs_square(wRl, wIl) == 0.) continue; // skip to next state
                    
                    // accumulate this layer amplitude onto this path amplitude
                    complex_multiply (path_wR, path_wI, acc_wR, acc_wI, wRl, wIl);
                    
                    // propagate via the very last layer towards final_state
                    layer_w(layer_last, l+1, next_state, final_state, wRll, wIll);
#if defined(DEBUG) && !defined(__FLOAT_MP__)
                    fprintf (stderr, "Layer l=%d -> %f + i %f\n", l+1, wRll, wIll);
#endif
                    // check whether this layer amplitude is zero
                    if (complex_abs_square(wRll, wIll) == 0.) continue; // skip to next state
                    
                    // build on the path amplitude
                    complex_multiply (path_wR, path_wI, path_wR, path_wI, wRll, wIll);
                    
                    // if this path amplitude is larger than 0 add it to final estimate
                    if (complex_abs_square(path_wR, path_wI) > 0.) {
                        sumR += path_wR;
                        sumI += path_wI;
                        
#if defined(DEBUG) && !defined(__FLOAT_MP__)
                         fprintf (stderr, "Path: %llu -> ", init_state);
                         for (int l=0 ; l<(L-2); l++) {
                         fprintf (stderr, "%llu -> ", ndxs[l]);
                           }
                         fprintf (stderr, "%llu -> %llu = ", next_state, final_state);
                         fprintf(stderr, "%f + i %f\n", path_wR, path_wI);
                         fprintf(stderr, "Current sum = %f + i %f\n\n", sumR, sumI);
#endif
                         
#ifdef NON_ZERO_PATHS
                        non_zero_paths++;
                        // output the path
                        T_Non_zero_path myPath;
                        myPath.path.push_back(init_state);
                        for (int l=0 ; l<L-2 ; l++) {
                            myPath.path.push_back(ndxs[l]);
                        }
                        myPath.path.push_back(next_state);
                        myPath.path.push_back(final_state);
                        myPath.path_ndx = path_hash(myPath.path, NBR_STATES);
                        myPath.Ar = path_wR;
                        myPath.Ai = path_wI;
                       
#if !defined(__FLOAT_MP__)
                        fprintf (stderr, "\t%llu : ", myPath.path_ndx);
                        for (auto& state: myPath.path) {
                            fprintf (stderr, "%llu -> ", state);
                        }
                        //for (int l=0 ; l<L-2 ; l++) fprintf (stderr, "%llu -> ", ndxs[l]);
                        //fprintf (stderr, "%llu -> %llu = ", next_state, final_state);
                        fprintf(stderr, "%f + i %f\n", myPath.Ar, myPath.Ai);
#endif
#endif
                    }
                }  // end paths
            } // end of processing the last 2 layers of gates
            
#ifdef DEBUG
            fprintf (stderr, "From here is control of for loops\n");
#endif
            // Control the remaining L-2 FOR loops' counters (in ndxs[])
            for (int l= L-3 ; l>=1 ; l--) {
                bool break_for_loop = false;
                
                // Adding 1 to this layer's state number will keep happening
                // if any of the associated layers (the layer befoire and the layer after)
                // have zero amplitude, since then the path also has 0 amplitude
                bool zero_power_transition = true;
                while (zero_power_transition and !break_for_loop and !stop) {
                    zero_power_transition = false;
                    ndxs[l]++;  // if execution gets here we do have to recompute this layer amplitude
                    if (ndxs[l] == NBR_STATES) {  // reset to 0 and proceeed to next for
                        if (l==1) {
                            stop=true;  // outer most loop: finish
                            break_for_loop = true;
                        }
                        else ndxs[l] = 0;
                    }
                    else break_for_loop = true;
                    
                    if (ndxs[l]!=0) { // do not do the |0> . It has been done above
                        myReal wR=1.0f, wI=0.0f;
                        TCircuitLayer *layer;
                        // Recompute amplitude through the NEXT gate layer  [l+1]
                        // current_state is ndxs[l]
                        // next_state is ndxs[l+1]  NOTE: only do this if l < L-3
                        if (l < L-3) {
                            // get gate layer l+1
                            layer = &c->layers[l+1];
                            // get the amplitude of this layer for the given current and next states
                            layer_w (layer, l+1, ndxs[l], ndxs[l+1], wR, wI);
                            
                            w_intermediateR[l+1] = wR;
                            w_intermediateI[l+1] = wI;

                            if (complex_abs_square(wR, wI) <= 0.) zero_power_transition = true;

                        }
                        // Recompute amplitude through the PREVIOUS gate layer  [l]
                        // current_state is ndxs[l-1]
                        // next_state is ndxs[l]
                            
                        // get gate layer l
                        layer = &c->layers[l];
                        // get the amplitude of this layer for the given current and next states
                        layer_w (layer, l, ndxs[l-1], ndxs[l], wR, wI);
                            
                        w_intermediateR[l] = wR;
                        w_intermediateI[l] = wI;
                        if (complex_abs_square(wR, wI) <= 0.) zero_power_transition = true;
                    }
                    
                } // while loop for consecutive increments of the same layer's state
                if (break_for_loop) break;  // no further propagation to lesser nested FOR loops
            }
#ifdef DEBUG
            fprintf (stderr, "END of control of for loops\n");
#endif
            
            
        }  // end of big loop implementing all the FOR loops
    } // outer most for loop (ndxs[0])
    
    estimateR = sumR;
    estimateI = sumI;
    
    delete [] w_intermediateI;
    delete [] w_intermediateR;
    delete [] ndxs;

    return true;
}



