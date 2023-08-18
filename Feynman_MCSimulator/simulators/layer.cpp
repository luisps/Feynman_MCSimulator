//
//  layer.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 18/08/2023.
//

#include "layer.hpp"
#include <random>

void layer_w (TCircuitLayer *layer, int l,
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
        const int qb_nbr = gate->fdata.qubit;
        const int curr_state_qb = qb_value(qb_nbr, current_state);
        const int next_state_qb = qb_value(qb_nbr, next_state);
        
        
        switch (gate->fdata.name) {
            case 11:            // rx
            case 12:            // ry
            case 13:            // rz
            case 14:            // p
                gate_g1p1_w (curr_state_qb, next_state_qb, gate->fdata.m, gatewR, gatewI);
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
        const int c_qb_nbr = gate->fdata.c_qubit;
        const int t_qb_nbr = gate->fdata.t_qubit;
        const int curr_state_qbs = qb_value(c_qb_nbr, current_state)*2+qb_value(t_qb_nbr, current_state);
        const int next_state_qbs = qb_value(c_qb_nbr, next_state)*2+qb_value(t_qb_nbr, next_state);
        
        
        switch (gate->fdata.name) {
            case 31:            // cp
                gate_g2p1_w (curr_state_qbs, next_state_qbs, gate->fdata.m, gatewR, gatewI);
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


float layer_sample (TCircuitLayer* layer, int l, unsigned long long current_state,
                    unsigned long long& next_state,
                    float& wR, float& wI,
                    std::default_random_engine e, std::uniform_real_distribution<float>& d) {
    
    float lwR = 1.f, lwI = 0.f, pdf=1.f;
    unsigned long long lnext_state=0ull;

    //if (print_debug) fprintf (stderr, "Evaluate layer %d with curr=%llu and next=%llu\n",l, current_state, next_state);
    // Iterate over all gates G1P0
    for (int g=0 ; g<layer->num_type_gates[0] ; g++) {
        float gatewR, gatewI, gatepdf;
        TGate1P0 *gate = &(((TGate1P0 *) layer->gates[0])[g]);
        const int qb_nbr = gate->qubit;
        const int curr_state_qb = qb_value(qb_nbr, current_state);
        //if (print_debug) fprintf (stderr, "\tcurr: qb_nbr = %d -> qb=%d\n",qb_nbr, curr_state_qb);
        const float rnd = d(e);  // generate random nbr in [0, 1[

        int next_state_qb = 0;

        switch (gate->name) {
            case 0:             // id
                gatepdf = gate_id_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
            case 1:             // h
                gatepdf = gate_h_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
            case 2:             // x
                gatepdf = gate_x_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
            case 3:             // y
                gatepdf = gate_y_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
            case 4:             // z
                gatepdf = gate_z_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
            case 5:             // s
                gatepdf = gate_s_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
            case 6:             // t
                gatepdf = gate_t_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
            default:             // id
                gatepdf = gate_id_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
        }
        //if (print_debug) fprintf (stderr, "\tnext: qb_nbr = %d -> qb=%d\n",qb_nbr, next_state_qb);
        //if (print_debug) fprintf (stderr, "Current Layer Product = %f + i %f\n", lwR, lwI);
        
        // verify early termination conditions
        if (complex_abs_square(gatewR, gatewI)<= 0.f || gatepdf <= 0.f) {
            return -1.f;
        }
        
        // Set the appropriate bit in lnext_state
        lnext_state = (next_state_qb==0 ? qb_reset (qb_nbr, lnext_state) : qb_set (qb_nbr, lnext_state));
        
        //if (print_debug) fprintf (stderr, "Layer %d, gate %d, name %d, qubit %d, in %d, out %d -> %f +i %f\n", l, g, gate->name, qb_nbr, curr_state_qb, next_state_qb, gatewR, gatewI);
        complex_multiply (lwR, lwI, lwR, lwI, gatewR, gatewI);
        pdf *= gatepdf;
        
        //if (print_debug) fprintf (stderr, "Layer Product = %f + i %f\n", lwR, lwI);
    }  // end all gates G1P0
    // Iterate over all gates G1P1
    for (int g=0 ; g<layer->num_type_gates[1] ; g++) {
        float gatewR, gatewI, gatepdf;
        TGate1P1 *gate = &(((TGate1P1 *) layer->gates[1])[g]);
        const int qb_nbr = gate->fdata.qubit;
        const int curr_state_qb = qb_value(qb_nbr, current_state);
        const float rnd = d(e);  // generate random nbr in [0, 1[

        int next_state_qb = 0;

        switch (gate->fdata.name) {
            case 11:            // rx
            case 12:            // ry
            case 13:            // rz
            case 14:            // p
                gatepdf = gate_g1p1_sample (curr_state_qb, rnd, next_state_qb, gate->fdata.m, gate->pdf, gatewR, gatewI);
                break;
            default:             // id
                gatepdf = gate_id_sample (curr_state_qb, rnd, next_state_qb, gatewR, gatewI);
                break;
        }
        // verify early termination conditions
        if (complex_abs_square(gatewR, gatewI)<= 0.f || gatepdf <= 0.f) {
            return -1.f;
        }
        
        // Set the appropriate bit in lnext_state
        lnext_state = (next_state_qb==0 ? qb_reset (qb_nbr, lnext_state) : qb_set (qb_nbr, lnext_state));

        complex_multiply (lwR, lwI, lwR, lwI, gatewR, gatewI);
        pdf *= gatepdf;

    }  // end all gates G1P1
    
    // Iterate over all gates G2P0
    for (int g=0 ; g<layer->num_type_gates[2] ; g++) {
        float gatewR, gatewI, gatepdf;
        TGate2P0 *gate = &(((TGate2P0 *) layer->gates[2])[g]);
        const int c_qb_nbr = gate->c_qubit;
        const int t_qb_nbr = gate->t_qubit;
        const int curr_state_qbs = qb_value(c_qb_nbr, current_state)*2+qb_value(t_qb_nbr, current_state);
        int next_state_qbs = 0;
        const float rnd = d(e);  // generate random nbr in [0, 1[
        
        switch (gate->name) {
            case 20:             // id2
                gatepdf = gate_id2_sample (curr_state_qbs, rnd, next_state_qbs, gatewR, gatewI);
                break;
            case 21:             // cx
                gatepdf = gate_cx_sample (curr_state_qbs, rnd, next_state_qbs, gatewR, gatewI);
                break;
            case 22:             // cz
                gatepdf = gate_cz_sample (curr_state_qbs, rnd, next_state_qbs, gatewR, gatewI);
                break;
            default:             // id2
                gatepdf = gate_id2_sample (curr_state_qbs, rnd, next_state_qbs, gatewR, gatewI);
                break;
        }
        // verify early termination conditions
        if (complex_abs_square(gatewR, gatewI)<= 0.f || gatepdf <= 0.f) {
            return -1.f;
        }
        
        // Set the appropriate bit in lnext_state
        lnext_state = ((next_state_qbs & 0x01) == 0 ? qb_reset (t_qb_nbr, lnext_state) : qb_set (t_qb_nbr, lnext_state));
        lnext_state = ((next_state_qbs & 0x02) == 0 ? qb_reset (c_qb_nbr, lnext_state) : qb_set (c_qb_nbr, lnext_state));

        complex_multiply (lwR, lwI, lwR, lwI, gatewR, gatewI);
        pdf *= gatepdf;
    } // end all gates G2P0
    
    // Iterate over all gates G2P1
    for (int g=0 ; g<layer->num_type_gates[3] ; g++) {
        float gatewR, gatewI, gatepdf;
        TGate2P1 *gate = &(((TGate2P1 *) layer->gates[3])[g]);
        const int c_qb_nbr = gate->fdata.c_qubit;
        const int t_qb_nbr = gate->fdata.t_qubit;
        const int curr_state_qbs = qb_value(c_qb_nbr, current_state)*2+qb_value(t_qb_nbr, current_state);
        const float rnd = d(e);  // generate random nbr in [0, 1[

        int next_state_qbs = 0;
        
        switch (gate->fdata.name) {
            case 31:            // cp
                gatepdf = gate_g2p1_sample (curr_state_qbs, rnd, next_state_qbs, gate->fdata.m, gate->pdf, gate->cdf, gatewR, gatewI);
                break;
            default:             // id
                gatepdf = gate_id2_sample (curr_state_qbs, rnd, next_state_qbs, gatewR, gatewI);
                break;
        }
        // verify early termination conditions
        if (complex_abs_square(gatewR, gatewI)<= 0.f || gatepdf <= 0.f) {
            return -1.f;
        }
        
        // Set the appropriate bit in lnext_state
        lnext_state = ((next_state_qbs & 0x01) == 0 ? qb_reset (t_qb_nbr, lnext_state) : qb_set (t_qb_nbr, lnext_state));
        lnext_state = ((next_state_qbs & 0x02) == 0 ? qb_reset (c_qb_nbr, lnext_state) : qb_set (c_qb_nbr, lnext_state));

        complex_multiply (lwR, lwI, lwR, lwI, gatewR, gatewI);
        pdf *= gatepdf;

    }
    // end all gates in the layer
    wR = lwR; wI = lwI;
    next_state = lnext_state;
    return pdf;
}
