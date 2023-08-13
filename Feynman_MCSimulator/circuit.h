//
//  circuit.h
//  Feynman_MCSimulator
//
//  Created by Luis Paulo Santos on 30/07/2023.
//

#ifndef circuit_h
#define circuit_h


typedef struct {
    int num_qubits;
    int num_layers;
} TCircuitSize;

typedef struct {
    int name;
    int qubit;
} TGate1P0;    // 1 qubit, 0 parameters gate datatype

typedef struct {
    int name;
    int qubit;
    float param;
    float m[2][2][2];  // note REAL and IMAG separated
} TGate1P1;    // 1 qubit, 1 parameter gate datatype

typedef struct {
    int name;
    int c_qubit, t_qubit;
    float param;
    float m[4][4][2];  // note REAL and IMAG separated
} TGate2P1;    // 2 qubits, 1 parameter gate datatype

typedef struct {
    int name;
    int c_qubit, t_qubit;
} TGate2P0;    // 2 qubits, 0 parameters gate datatype

typedef struct {
    int num_gates;
    int num_type_gates[4];  // G1P0,G1P1,G2P0,G2P1
    void *gates[4];         // G1P0,G1P1,G2P0,G2P1
} TCircuitLayer;


typedef struct {
    TCircuitSize *size;
    TCircuitLayer *layers;
} TCircuit;

TCircuit * read_circuit (const char *);

void print_circuit (TCircuit *);

#endif /* circuit_h */
