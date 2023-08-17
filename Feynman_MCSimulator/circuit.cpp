//
//  circuit.cpp
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 30/07/2023.
//

#include <stdio.h>
#include <stdlib.h>

#include "circuit.h"
#include "complex.h"


TCircuit * read_circuit (const char *fileName) {
    FILE *f=NULL;
    TCircuit *circuit=NULL;
    TCircuitSize *csize=NULL;
    TCircuitLayer *clayers=NULL;
    int fread_ret;
    
    
    f = fopen(fileName, "rb");
    if (f==NULL) {
        fprintf (stderr, "read_circuit() could not open file %s\n", fileName);
        return NULL;
    }
    
    circuit = (TCircuit *) malloc(sizeof(TCircuit));
    if (circuit==NULL) {
        fprintf (stderr, "read_circuit() could not malloc mem for circuit\n");
        return NULL;
    }
    csize = (TCircuitSize *) malloc(sizeof(TCircuitSize));
    if (csize == NULL) {
        fprintf (stderr, "read_circuit() could not malloc mem for circuit\n");
        return NULL;
    }

    // read circuit height and width
    fread_ret = (int) fread((void *)csize, sizeof(TCircuitSize), 1, f);
    if ( fread_ret != 1) {
        fprintf (stderr, "read_circuit() could not read TCircuitSize structure\n");
        return NULL;
    }
    circuit->size = csize;
    
    // alloc mem for layers
    clayers = (TCircuitLayer *) malloc(csize->num_layers * sizeof(TCircuitLayer));
    if (clayers == NULL) {
        fprintf (stderr, "read_circuit() could not malloc mem for layers\n");
        return NULL;
    }
    circuit->layers = clayers;
    
    // for each layer
    for (int l=0 ; l<csize->num_layers ; l++) {
        TCircuitLayer *layer = &clayers[l];
        struct ff {int G1P0, G1P1, G2P0, G2P1;} data;
        // read layer's number of different gates
        fread_ret = (int) fread((void *)&data, sizeof(data), 1, f);
        if ( fread_ret != 1) {
            fprintf (stderr, "read_circuit() could not read layer %d gate arrangement\n", l);
            return NULL;
        }
        
        // Number of gates for this layer
        layer->num_type_gates[0] = data.G1P0;
        layer->num_type_gates[1] = data.G1P1;
        layer->num_type_gates[2] = data.G2P0;
        layer->num_type_gates[3] = data.G2P1;
        layer->num_gates = data.G1P0 + data.G1P1 + data.G2P0 + data.G2P1;
        
        // read in G1P0
        if (layer->num_type_gates[0] > 0) {
            // alloc mem for layers[l].num_type_gates[0] gates
            layer->gates[0] = (void *)malloc(layer->num_type_gates[0] * sizeof(TGate1P0));
            if (layer->gates[0] == NULL) {
                fprintf (stderr, "read_circuit() could not malloc mem for G1P0 gates\n");
                return NULL;
            }
            TGate1P0 * g1p0_ptr = (TGate1P0 *) layer->gates[0];
            for (int g=0 ; g<layer->num_type_gates[0] ; g++) {
                // read next G1P0 gate
                fread_ret = (int) fread((void *)g1p0_ptr, sizeof(TGate1P0), 1, f);
                if ( fread_ret != 1) {
                    fprintf (stderr, "read_circuit() could not read layer %d gate G1P0 nbr %d\n", l, g);
                    return NULL;
                }
                //fprintf(stderr, "Read in gate name %d and qubit %d\n", g1p0_ptr->name, g1p0_ptr->qubit);
                g1p0_ptr++;
                
            }
        }

        // read in G1P1
        
        if (layer->num_type_gates[1] > 0) {
            // alloc mem for layers[l].num_type_gates[1] gates
            layer->gates[1] = (void *)malloc(layer->num_type_gates[1] * sizeof(TGate1P1));
            if (layer->gates[1] == NULL) {
                fprintf (stderr, "read_circuit() could not malloc mem for G1P1 gates\n");
                return NULL;
            }
            TGate1P1 * g1p1_ptr = (TGate1P1 *) layer->gates[1];
            for (int g=0 ; g<layer->num_type_gates[1] ; g++) {
                // read next G1P1 gate: only part of the data is on file
                fread_ret = (int) fread((void *)&g1p1_ptr->fdata, sizeof(TGate1P1_FDATA), 1, f);
                if ( fread_ret != 1) {
                    fprintf (stderr, "read_circuit() could not read layer %d gate G1P1 nbr %d\n", l, g);
                    return NULL;
                }
                // compute the pdf data
                for (int r=0 ; r<2 ; r++) {
                    for (int c=0 ; c<2 ; c++) {
                        g1p1_ptr->pdf[r][c] = complex_abs_square(g1p1_ptr->fdata.m[r][c][0],g1p1_ptr->fdata.m[r][c][1]);
                    }
                }
                //fprintf(stderr, "Read in gate name %d and qubit %d\n", g1p0_ptr->name, g1p0_ptr->qubit);
                g1p1_ptr++;
            }
        }

        // read in G2P0
        
        if (layer->num_type_gates[2] > 0) {
            // alloc mem for layers[l].num_type_gates[0] gates
            layer->gates[2] = (void *)malloc(layer->num_type_gates[2] * sizeof(TGate2P0));
            if (layer->gates[2] == NULL) {
                fprintf (stderr, "read_circuit() could not malloc mem for G2P0 gates\n");
                return NULL;
            }
            TGate2P0 * g2p0_ptr = (TGate2P0 *) layer->gates[2];
            for (int g=0 ; g<layer->num_type_gates[2] ; g++) {
                // read next G2P0 gate
                fread_ret = (int) fread((void *)g2p0_ptr, sizeof(TGate2P0), 1, f);
                if ( fread_ret != 1) {
                    fprintf (stderr, "read_circuit() could not read layer %d gate G2P0 nbr %d\n", l, g);
                    return NULL;
                }
                g2p0_ptr++;
            }
        }

        // read in G2P1
        
        if (layer->num_type_gates[3] > 0) {
        // alloc mem for layers[l].num_type_gates[3] gates
        layer->gates[3] = (void *)malloc(layer->num_type_gates[3] * sizeof(TGate2P1));
        if (layer->gates[3] == NULL) {
            fprintf (stderr, "read_circuit() could not malloc mem for G2P1 gates\n");
            return NULL;
        }
        TGate2P1 * g2p1_ptr = (TGate2P1 *) layer->gates[3];
            for (int g=0 ; g<layer->num_type_gates[3] ; g++) {
                // read next G2P1 gate: only part of the data is on file
                fread_ret = (int) fread((void *)&g2p1_ptr->fdata, sizeof(TGate2P1_FDATA), 1, f);
                if ( fread_ret != 1) {
                    fprintf (stderr, "read_circuit() could not read layer %d gate G2P1 nbr %d\n", l, g);
                    return NULL;
                }
                //fprintf(stderr, "Read in gate name %d and qubit %d\n", g1p0_ptr->name, g1p0_ptr->qubit);
                // compute the pdf and cdf data
                // NOTE: the cdf is transposed such that sampling an output for a given input uses a single row
                //       transposition occurs HERE in read_circuit()

                float acdf[4] = {0.f, 0.f, 0.f, 0.f};
                for (int r=0 ; r<4 ; r++) {
                    for (int c=0 ; c<4 ; c++) {
                        const float abs_sq = complex_abs_square(g2p1_ptr->fdata.m[r][c][0],g2p1_ptr->fdata.m[r][c][1]);
                        g2p1_ptr->pdf[r][c] = abs_sq;
                        acdf[c] += abs_sq;
                        // TRANSPOSITION
                        g2p1_ptr->cdf[c][r] = acdf[c];
                    }
                }
                g2p1_ptr++;
            }
        }

    }
    
    
    return circuit;
}

void print_circuit_stats (TCircuit *c) {
    fprintf (stdout, "Qubits: %d\n", c->size->num_qubits);
    fprintf (stdout, "Layers: %d\n", c->size->num_layers);
    fprintf(stdout, "\n");
}

void print_circuit (TCircuit *c) {

    print_circuit_stats(c);
    // for each layer
    for (int l=0 ; l < c->size->num_layers ; l++) {
        TCircuitLayer *layer= &c->layers[l];
        fprintf(stdout, "Layer %d\n", l);
        fprintf(stdout, "\t%d gates : G1P0 %d; G1P1 %d; G2P0 %d; G2P1 %d\n", layer->num_gates, layer->num_type_gates[0], layer->num_type_gates[1], layer->num_type_gates[2], layer->num_type_gates[3]);

        // for each gate type
        // G1P0  G1P1  G2P0  G2P1
        
        if (layer->num_type_gates[0] > 0) {  // G1P0
            fprintf(stdout, "\t\tG1P0\n");
            TGate1P0 * g1p0_ptr = (TGate1P0 *) layer->gates[0];
            for (int g=0 ; g< layer->num_type_gates[0] ; g++, g1p0_ptr++) {
                fprintf(stdout, "\t\tGate name %d; qubit %d\n", g1p0_ptr->name, g1p0_ptr->qubit);
            }
        }
        if (layer->num_type_gates[1] > 0) {  // G1P1
            fprintf(stdout, "\t\tG1P1\n");
            TGate1P1 * g1p1_ptr = (TGate1P1 *) layer->gates[1];
            for (int g=0 ; g< layer->num_type_gates[1] ; g++, g1p1_ptr++) {
                fprintf(stdout, "\t\tGate name %d; qubit %d; param %f\n", g1p1_ptr->fdata.name, g1p1_ptr->fdata.qubit, g1p1_ptr->fdata.param);
                for (int r=0 ; r<2 ; r++)
                    for (int c=0 ; c<2 ; c++)
                        fprintf(stdout, "\t\t\t[%d][%d] = %f + i %f\n", r,c,g1p1_ptr->fdata.m[r][c][0],g1p1_ptr->fdata.m[r][c][1]);
            }
        }
        if (layer->num_type_gates[2] > 0) {  // G2P0
            fprintf(stdout, "\t\tG2P0\n");
            TGate2P0 * g2p0_ptr = (TGate2P0 *) layer->gates[2];
            for (int g=0 ; g< layer->num_type_gates[2] ; g++, g2p0_ptr++) {
                fprintf(stdout, "\t\tGate name %d; C qubit %d; T qubit %d\n", g2p0_ptr->name, g2p0_ptr->c_qubit, g2p0_ptr->t_qubit);
            }
        }
        if (layer->num_type_gates[3] > 0) {  // G2P1
            fprintf(stdout, "\t\tG2P1\n");
            TGate2P1 * g2p1_ptr = (TGate2P1 *) layer->gates[3];
            for (int g=0 ; g< layer->num_type_gates[3] ; g++, g2p1_ptr++) {
                fprintf(stdout, "\t\tGate name %d; C qubit %d; T qubit %d; param %f\n", g2p1_ptr->fdata.name, g2p1_ptr->fdata.c_qubit, g2p1_ptr->fdata.t_qubit, g2p1_ptr->fdata.param);
                for (int r=0 ; r<4 ; r++)
                    for (int c=0 ; c<4 ; c++)
                        fprintf(stdout, "\t\t\t[%d][%d] = %f + i %f\n", r,c,g2p1_ptr->fdata.m[r][c][0],g2p1_ptr->fdata.m[r][c][1]);
            }
        }
    }
    fprintf(stdout, "\n");

}
