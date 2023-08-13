//
//  all_paths.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 31/07/2023.
//

#include <stdlib.h>

#include "all_paths.hpp"
#include "gates.h"
#include "complex.h"

static bool all_paths_NOVEC_NOOMP (TCircuit *, unsigned long long init_state,
                unsigned long long final_state,
                float &estimateR, float &estimateI);

static void layer_w (TCircuitLayer *, int ,
                     unsigned long long , unsigned long long , float &, float &, bool print_debug=false);

bool all_paths (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state,
                float &estimateR, float &estimateI) {
    
    bool ret;
    
    ret = all_paths_NOVEC_NOOMP (c, init_state, final_state, estimateR, estimateI);
    
    return ret;
}


static bool all_paths_NOVEC_NOOMP (TCircuit *c, unsigned long long init_state,
                unsigned long long final_state,
                                   float &estimateR, float &estimateI) {
    
    float sumR=0., sumI=0.;
    const int L = c->size->num_layers;  // number of gate layers
    /* the number of intermediate states layers is L-1.
       This is the number of for loops we have to realize
       We will do L-2 for loops in a single aggregated while loop
       The most nested for loop (corresponding to L-1 (just before trransitioning
       to the final state) is still implememnted with a for loop */
    
    /* Note that the loop counter for each intermediate state layer encodes
     the state itself. */
    
    unsigned long long *ndxs = (unsigned long long *) malloc ((L-2)* sizeof(unsigned long long ));
    
    // compute NBR_STATES : the number of states to iterate in each intermediate state layer
    const unsigned long long NBR_STATES = 1 << c->size->num_qubits;
    
    // initialize all our L-2 loop counters
    for (int l=0 ; l< (L-2) ; l++) ndxs[l] = 0ull;
    
    /* We will enter the multiple variable number of for loops below */
    /* We must first make all initializations corresponding to all ndxs==0 */
    /* Note that this corresponds to computing the amplitude of the path
       transitioning from the initial state, through all intermediate states = |0>
       all way to the final state.
       However, since the final transition to the final state is supported by an
       explicit for loop we don't pre compute the very last one */
    
    // store the intermediate amplitudes associated with each transition between layers
    float *w_intermediateR,  *w_intermediateI;
    w_intermediateR = (float *)malloc((L-2)*sizeof(float));
    w_intermediateI = (float *)malloc((L-2)*sizeof(float));
    unsigned long long current_state = init_state;
    
    for (int l=0 ; l < (L-2) ; l++) {
        float wR, wI;
        
        // get gate layer l
        TCircuitLayer *layer = &c->layers[l];
        
        // get the amplitude of this layer for the given current and next states
        layer_w (layer, l, current_state, ndxs[l], wR, wI);
        
        //fprintf (stderr, "Layer l=%d -> %f + i %f\n", l, wR, wI);
        
        w_intermediateR[l] = wR;
        w_intermediateI[l] = wI;

        // update current state
        current_state = ndxs [l];
    }  // end iterating over layers 0 .. L-3 (that is L-2 layers)

    /* the code corresponding to the first iteration of outer nested L-2 for loops
       has been executed
        Now enter the loops , starting with the inner most for loop, i.e.,
        from current_state ndxs[L-3] ->  all s in last intermediate state and than to final_state  \=>
     */
    
    bool stop = false;
    int non_zero_paths =0;
    while (!stop) {
        float acc_wR, acc_wI;
        acc_wR = 1.f; acc_wI=0.f;
        
        // check whether the amplitude from previous layers is 0
        for (int l=0 ; l<(L-2) ; l++) {
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
            
            current_state = (L >= 3 ? ndxs[L-3] : init_state);
            for (unsigned long long next_state=0ull; next_state<NBR_STATES ; next_state++) {
                float path_wR, path_wI;
                float wRl, wIl, wRll, wIll;
                
                layer_w(layer_2L, l, current_state, next_state, wRl, wIl);
                //fprintf (stderr, "Layer l=%d -> %f + i %f\n", l, wRl, wIl);
                // check whether this layer amplitude is zero
                if (complex_abs_square(wRl, wIl) == 0.) continue; // skip to next state
                
                // accumulate this layer amplitude onto this path amplitude
                complex_multiply (path_wR, path_wI, acc_wR, acc_wI, wRl, wIl);
                
                // propagate via the very last layer towards final_state
                layer_w(layer_last, l+1, next_state, final_state, wRll, wIll);
                //fprintf (stderr, "Layer l=%d -> %f + i %f\n", l+1, wRll, wIll);
                // check whether this layer amplitude is zero
                if (complex_abs_square(wRll, wIll) == 0.) continue; // skip to next state
                
                // build on the path amplitude
                complex_multiply (path_wR, path_wI, path_wR, path_wI, wRll, wIll);
                
                // if this path amplitude is larger than 0 add it to final estimate
                if (complex_abs_square(path_wR, path_wI) > 0.) {
                    sumR += path_wR;
                    sumI += path_wI;
                    
                    /*
                    fprintf (stderr, "Path: %llu -> ", init_state);
                    for (int l=0 ; l<(L-2); l++) {
                        fprintf (stderr, "%llu -> ", ndxs[l]);
                    }
                    fprintf (stderr, "%llu -> %llu = ", next_state, final_state);
                    fprintf(stderr, "%f + i %f\n", path_wR, path_wI);
                    fprintf(stderr, "Current sum = %f + i %f\n\n", sumR, sumI);
                     */
                    non_zero_paths++;
                    
                }
            }  // end paths
        } // end of processing the last 2 layers of gates
        
        //fprintf (stderr, "From here is control of for loops\n");
        // Control the remaining L-2 FOR loops' counters (in ndxs[])
        for (int l= L-3 ; l>=0 ; l--) {
            bool break_loop = false;
            ndxs[l]++;  // if execution gets here we do have to recompute this layer amplitude
            if (ndxs[l] == NBR_STATES) {  // reset to 0 and proceeed to next for
                if (l==0) {
                    stop=true;  // outer most loop: finish
                    break;
                }
                else ndxs[l] = 0;
            }
            else break_loop = true;
            
            float wR, wI;
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
            }
            // Recompute amplitude through the PREVIOUS gate layer  [l]
            // current_state is ndxs[l-1] or init_state
            // next_state is ndxs[l]

            // get gate layer l
            layer = &c->layers[l];
            // get the amplitude of this layer for the given current and next states
            layer_w (layer, l, (l==0 ? init_state : ndxs[l-1]), ndxs[l], wR, wI);
                
            w_intermediateR[l] = wR;
            w_intermediateI[l] = wI;

            if (break_loop) break;  // no further propagation to lesser nested FOR loops
        }
        //fprintf (stderr, "END of control of for loops\n");

        
    }  // end of big loop implementing all the FOR loops
    
    estimateR = sumR;
    estimateI = sumI;
    
    fprintf (stderr, "There were %d non zero paths.\n", non_zero_paths);
    
    return true;
}


static void layer_w (TCircuitLayer *layer, int l,
                     unsigned long long current_state, unsigned long long next_state, float &wR, float &wI, bool print_debug) {
    float lwR = 1.f, lwI = 0.f;

    //if (print_debug) fprintf (stderr, "Evaluate layer %d with curr=%llu and next=%llu\n",l, current_state, next_state);
    // Iterate over all gates G1P0
    for (int g=0 ; g<layer->num_type_gates[0] ; g++) {
        float gatewR, gatewI;
        TGate1P0 *gate = &(((TGate1P0 *) layer->gates[0])[g]);
        const int qb_nbr = gate->qubit;
        const int curr_state_qb = qb_value(qb_nbr, current_state);
        //if (print_debug) fprintf (stderr, "\tcurr: qb_nbr = %d -> qb=%d\n",qb_nbr, curr_state_qb);

        const int next_state_qb = qb_value(qb_nbr, next_state);
        //if (print_debug) fprintf (stderr, "\tnext: qb_nbr = %d -> qb=%d\n",qb_nbr, next_state_qb);
        //if (print_debug) fprintf (stderr, "Current Layer Product = %f + i %f\n", lwR, lwI);

        switch (gate->name) {
            case 0:             // id
                gate_id_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
            case 1:             // h
                gate_h_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
            case 2:             // x
                gate_x_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
            case 3:             // y
                gate_y_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
            case 4:             // z
                gate_z_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
            case 5:             // s
                gate_s_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
            case 6:             // t
                gate_t_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
            default:             // id
                gate_id_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
        }
        //if (print_debug) fprintf (stderr, "Layer %d, gate %d, name %d, qubit %d, in %d, out %d -> %f +i %f\n", l, g, gate->name, qb_nbr, curr_state_qb, next_state_qb, gatewR, gatewI);
        complex_multiply (lwR, lwI, lwR, lwI, gatewR, gatewI);
        //if (print_debug) fprintf (stderr, "Layer Product = %f + i %f\n", lwR, lwI);
    }  // end all gates G1P0
    // Iterate over all gates G1P1
    for (int g=0 ; g<layer->num_type_gates[1] ; g++) {
        float gatewR, gatewI;
        TGate1P1 *gate = &(((TGate1P1 *) layer->gates[1])[g]);
        const int qb_nbr = gate->qubit;
        const int curr_state_qb = qb_value(qb_nbr, current_state);
        const int next_state_qb = qb_value(qb_nbr, next_state);
        
        
        switch (gate->name) {
            case 11:            // rx
            case 12:            // ry
            case 13:            // rz
            case 14:            // p
                gate_g1p1_w (curr_state_qb, next_state_qb, gate->m, gatewR, gatewI);
                break;
            default:             // id
                gate_id_w (curr_state_qb, next_state_qb, gatewR, gatewI);
                break;
        }
        complex_multiply (lwR, lwI, lwR, lwI, gatewR, gatewI);
    }  // end all gates G1P1
    
    // Iterate over all gates G2P0
    for (int g=0 ; g<layer->num_type_gates[2] ; g++) {
        float gatewR, gatewI;
        TGate2P0 *gate = &(((TGate2P0 *) layer->gates[2])[g]);
        const int c_qb_nbr = gate->c_qubit;
        const int t_qb_nbr = gate->t_qubit;
        const int curr_state_qbs = qb_value(c_qb_nbr, current_state)*2+qb_value(t_qb_nbr, current_state);
        const int next_state_qbs = qb_value(c_qb_nbr, next_state)*2+qb_value(t_qb_nbr, next_state);
        switch (gate->name) {
            case 20:             // id2
                gate_id2_w (curr_state_qbs, next_state_qbs, gatewR, gatewI);
                break;
            case 21:             // cx
                gate_cx_w (curr_state_qbs, next_state_qbs, gatewR, gatewI);
                break;
            case 22:             // cz
                gate_cz_w (curr_state_qbs, next_state_qbs, gatewR, gatewI);
                break;
            default:             // id2
                gate_id2_w (curr_state_qbs, next_state_qbs, gatewR, gatewI);
                break;
        }
        complex_multiply (lwR, lwI, lwR, lwI, gatewR, gatewI);
    } // end all gates G2P0
    
    // Iterate over all gates G2P1
    for (int g=0 ; g<layer->num_type_gates[3] ; g++) {
        float gatewR, gatewI;
        TGate2P1 *gate = &(((TGate2P1 *) layer->gates[3])[g]);
        const int c_qb_nbr = gate->c_qubit;
        const int t_qb_nbr = gate->t_qubit;
        const int curr_state_qbs = qb_value(c_qb_nbr, current_state)*2+qb_value(t_qb_nbr, current_state);
        const int next_state_qbs = qb_value(c_qb_nbr, next_state)*2+qb_value(t_qb_nbr, next_state);
        
        
        switch (gate->name) {
            case 31:            // cp
                gate_g2p1_w (curr_state_qbs, next_state_qbs, gate->m, gatewR, gatewI);
                break;
            default:             // id
                gate_id2_w (curr_state_qbs, next_state_qbs, gatewR, gatewI);
                break;
        }
        complex_multiply (lwR, lwI, lwR, lwI, gatewR, gatewI);
    }
    // end all gates in the layer
    wR = lwR; wI = lwI;
}
