# generates random circuits as described in Peres, Filipa 2023 (sec 4.2)
# https://arxiv.org/pdf/2203.01789.pdf

import qiskit
import random


### To remove pairs of successive Hadamards from a qiskit circuit

def get_gate_qubits (gate, qc):
    g_qubits = []
    for qb in gate[1]:
        g_qubits.append(qc.find_bit(qb).index)
    return g_qubits
    
def remove_H_pairs (qc):
    ndx_to_remove = []
    n = qc.num_qubits;
    # lists with H gate info per qubit
    lastH = [False] * n
    ndx_lastH = [0] * n
    
    # iterate over gates
    for ndx, gate in enumerate(qc.data):
        name = gate[0].name
        g_qubits = get_gate_qubits (gate, qc)
        if name != 'h':  # not an Hadamard
            for qb in g_qubits:
                lastH[qb] = False
        else:            # Hadamard
            qb = g_qubits[0]    # single qubit gate
            if not lastH[qb]:   # first Hadamard of a potential pair
                lastH[qb] = True
                ndx_lastH[qb] = ndx
            else:               # second Hadamard in a pair
                ndx_to_remove.append (ndx)
                ndx_to_remove.append (ndx_lastH[qb])
                lastH[qb] = False
    
    #print (ndx_to_remove)
    for i in sorted(ndx_to_remove, reverse=True):
        del qc.data[i]
    return qc
### END:: To remove pairs of successive Hadamards from a qiskit circuit

###################
#
# Exported function
#
#. Note: the number of qubits must be even

def ran_gen (n, nCycles, measure = True, remove_pairs_H=False):
        
    qr = qiskit.QuantumRegister(n, 'qr')
    if measure:
        cr = qiskit.ClassicalRegister(n, 'cr')
        qc = qiskit.QuantumCircuit(qr, cr, name='U')
    else:
        qc = qiskit.QuantumCircuit(qr, name='U')
        
    # Initial layer of Hadamards
    for i in range(n):
        qc.h(i)
        
    # nCycles (layers)
    # each layer connects n/4 pairs of qubits using CZ (H - CX - H) - these are n/2 qubits
    aHalf = n // 2
    aHalf += 0 if aHalf%2==0 else 1 #make it even
    for layer in range (nCycles):
        toEntangle = random.sample (range(n), aHalf)
        for qb_ndx in range(0, aHalf, 2):
            qb_c = toEntangle[qb_ndx]
            qb_t = toEntangle[qb_ndx+1]
            qc.h(qb_t)
            qc.cx (qb_c, qb_t)
            qc.h(qb_t)

        # apply H,S or T on non connected qubits
        nonConnected = [qb for qb in range(n)]
        # remove the connected ones
        for i in sorted(toEntangle, reverse=True):
            del nonConnected[i]
            
        # Gates (see below) will be randomly selected among the following
        # 0 - 'h', 1 - 's', 2 - 't'
        for qb in nonConnected:
            gate = random.randrange(3)
            if gate == 0:
                    qc.h(qb)
                    
            elif gate == 1:
                    qc.s(qb)
                    
            elif gate == 2:
                    qc.t(qb)


    # END iterating nCycles (layers)
    if measure:
        for i in range(n):
            qc.measure(i, i)

    if remove_pairs_H:
        qc = remove_H_pairs (qc)
        
    return qc