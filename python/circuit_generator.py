from qiskit import QuantumCircuit
from qiskit.circuit.library import IQP, EfficientSU2
from qiskit.circuit.parametervector import ParameterVector, ParameterVectorElement
from qiskit.circuit import Parameter
from qiskit import quantum_info as qiskit_info

from HSA import HSA_gen
from RandomFilipa2023 import ran_gen

import random
from math import pi

##### Circuits description
#
# 1 - qbs=3, depth =3; [HHT][THS][HIH]4/64 paths !=0; constructive interference
# 2 - qbs=4, depth =3; [HITH][THSH][HHHH]16/256 paths !=0; but 8/16 paths self anhilate
# 3 - qbs=3, depth =3; [HHH][HHH][HHH]64/64 paths !=0; but almost half self anhilate due to destructive interference
# 4 - qbs=2, depth =2; [HH][HH] 4/4 paths !=0; only constructive interference
# 5 - qbs=2, depth =4; [HH][HH][SS] [HH] 
# 6 - qbs=2, depth =4; [HH][HH][TT][HH] 16/64 paths !=0, but 12/16 paths self anhilate
# 7 - qbs=3, depth =4; [HHHZ][HHHZ][HHHZ]64/512 paths !=0; but almost half self anhilate due to destructive interference
# 8 - qbs=4, depth =4; [HHHH][HHHH][HHHH][HHHH] 4096/4096 paths !=0; but almost half self anhilate 
# 9 - 

### To remove pairs of successive Hadamards from a qiskit circuit

def _get_gate_qubits (gate, qc):
    g_qubits = []
    for qb in gate[1]:
        g_qubits.append(qc.find_bit(qb).index)
    return g_qubits
  
def _remove_H_pairs (qc):
    ndx_to_remove = []
    n = qc.num_qubits;
    # lists with H gate info per qubit
    lastH = [False] * n
    ndx_lastH = [0] * n
    
    # iterate over gates
    for ndx, gate in enumerate(qc.data):
        name = gate[0].name
        g_qubits = _get_gate_qubits (gate, qc)
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

# converts a QuantumCircuit to our internal layered representation
def QCircuit_to_layers (qc):
    num_qubits = qc.num_qubits
    qc_data = qc.data

    #print ('Circuit has {0} qubits'.format(num_qubits))
    #print ()
    #print (qc_data)


    curr_layer = [0] * num_qubits  # current layer for each qubit
    layers = []
    for g in qc_data:  # each gate in the circuit
        num_layers = len(layers)
        #print ('Circuit has {0} at this stage.'.format(num_layers)) 
        instruction = g[0]
        n_qubits = instruction.num_qubits
        gate = instruction.name
        params = instruction.params
        #derreferenced_params = []
        #for param in params:
        #    if isinstance(param, ParameterVectorElement):
        #        #print (type(param), flush=True)
        #        #print (param)
        #        #print (params_data)
        #        #print (param.index)
        #        derreferenced_params.append(params_data[param])
        #        print (derreferenced_params[-1])
        #    else:
        #        derreferenced_params.append(param)
        # gates which depend on parameters have their unitary precomputed
        #print(params)
        #print (g)
        pdfs = None
        unitary = None
        if gate in ['rx','ry','rz','p','cp']:
            # Build a circuit to evaluate the parameterized unitary
            thisGateqc = QuantumCircuit(n_qubits)
            if gate=='rx':  
                thisGateqc.rx(params[0],0)
                op = qiskit_info.Operator (thisGateqc)
                unitary = op.data
                pdfs = []
                for row in [0,1]:
                    pdfs.append([])
                    for col in [0,1]:
                        unit = unitary[row][col]
                        prob_c = unit*unit.conjugate()
                        prob = prob_c.real
                        pdfs[-1].append(prob)
            elif gate=='ry':  
                thisGateqc.ry(params[0],0)
                op = qiskit_info.Operator (thisGateqc)
                unitary = op.data
                pdfs = []
                for row in [0,1]:
                    pdfs.append([])
                    for col in [0,1]:
                        unit = unitary[row][col]
                        prob_c = unit*unit.conjugate()
                        prob = prob_c.real
                        pdfs[-1].append(prob)
            elif gate=='rz':  
                thisGateqc.rz(params[0],0)
                op = qiskit_info.Operator (thisGateqc)
                unitary = op.data
                pdfs = []
                for row in [0,1]:
                    pdfs.append([])
                    for col in [0,1]:
                        unit = unitary[row][col]
                        prob_c = unit*unit.conjugate()
                        prob = prob_c.real
                        pdfs[-1].append(prob)
            elif gate=='p':  
                thisGateqc.p(params[0],0)
                op = qiskit_info.Operator (thisGateqc)
                unitary = op.data
            elif gate=='cp':  
                thisGateqc.cp(params[0],0,1)
                op = qiskit_info.Operator (thisGateqc)
                unitary = op.data
            #print ('Parameterized {0} gate with unitary: {1}'.format(gate,unitary))
        #print ('Gate {0} acts on {1} qubits: '.format(gate, n_qubits), end='')
        #print (g)
        layer_to_use = 0
        gate_qubits = []
        for qb in g[1]:
            qb_index = qb.index
            gate_qubits.append(qb_index)
            #print ('{0} '.format(qb_index),end='')
            if curr_layer[qb_index] >= layer_to_use:
                layer_to_use = curr_layer[qb_index]
            #print ('qubit {0} requires layer {1}'.format(qb_index, layer_to_use))
        if layer_to_use >= num_layers:  # add new layer
            layers.append([None]*num_qubits)
        # use this layer
        layer = layers[layer_to_use]
        for qb in gate_qubits:
            #layer[qb] = ((gate, gate_qubits, derreferenced_params, unitary, pdfs))
            layer[qb] = ((gate, gate_qubits, params, unitary, pdfs))
            curr_layer[qb]=layer_to_use+1 
        #print ('curr layers after this gate is ', curr_layer)
    # fill in None with 'id'
    for layer in layers:
        for qb in range(num_qubits):
            if layer[qb] is None:
                layer[qb] = (('id',[qb],[],[],None))
    return layers, len(layers)

def circuit1 ():
    qc = QuantumCircuit(3)
    #Layer 0
    qc.h(0) 
    #qc.id(1) 
    qc.t(2)
    #Layer 1
    qc.t(0)
    qc.h(1)
    qc.s(2)
    #Layer 2
    qc.cnot(1,0)
    #qc.h(0)
    #qc.h(1)
    qc.h(2)
    return qc    

def circuit2 ():
    qc = QuantumCircuit(4)
    #Layer 0
    qc.h(0) 
    qc.id(1) 
    qc.t(2)
    qc.h(3)
    #Layer 1
    qc.t(0)
    qc.h(1)
    qc.s(2)
    qc.h(3)
    #Layer 2
    qc.h(0)
    qc.h(1)
    qc.h(2)
    qc.h(3)
    #Layer 3
    qc.h(0) 
    qc.h(1)
    qc.h(2)
    qc.id(3) 
    return qc   

def circuit3 (num_layers, num_qubits):
    qc = QuantumCircuit(num_qubits)
    for l in range(num_layers):
        for qb in range(num_qubits):
            qc.h(qb) 
    return qc   

def circuit6 (num_qubits):
    qc = QuantumCircuit(num_qubits)
    for l in range(4):
        for qb in range(num_qubits):
            if l==2:
                qc.t(qb)
            else:
                qc.h(qb) 
    return qc   

def circuit5 (num_qubits):
    qc = QuantumCircuit(num_qubits)
    for l in range(3):
        for qb in range(num_qubits):
            if l==1:
                qc.t(qb)
            else:
                qc.h(qb) 
    return qc   


def get_circuit(ID, *args):
    layers = []

    if ID == 1:
        qc = circuit1()
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==2:
        qc = circuit2()
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers

    elif ID==3:
        num_qubits = args[0]
        num_layers = num_qubits
        qc = circuit3(num_layers, num_qubits)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
            
    elif ID==4:
        num_qubits = 2
        num_layers = 2
        qc = circuit3(num_layers, num_qubits)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers

    elif ID==5:
        qc = circuit5(4)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==51:
        qc = circuit5(1)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==52:
        qc = circuit5(2)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers

    elif ID==6:
        qc = circuit6(4)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==61:
        qc = circuit6(1)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==62:
        qc = circuit6(2)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers

    elif ID==7:
        num_qubits = 3
        num_layers = 4
        qc = circuit3(num_layers, num_qubits)
        for qb in range(num_qubits):
            qc.z(qb)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers

    elif ID==8:
        num_qubits = 4
        num_layers = 4
        qc = circuit3(num_layers, num_qubits)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers

            
    elif ID==9:
        A = [[6, 5, 3, 2], [5, 4, 5, 3], [3, 5, 2, 1], [2, 3, 1, 3]]
        qc = IQP(A)
        qc = qc.decompose()
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers

    elif ID==10:
        num_qubits = 2
        qc = circuit3(1, num_qubits)
        qc.cz(0,1)             
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers

    elif ID==11:
        num_qubits = 3
        qc = QuantumCircuit(num_qubits)
        qc.rx(pi/4,0)             
        qc.rx(pi/2,1)  
        qc.rx(pi,2)             
        qc.cx (2,1)      
        qc.cx (0,1)
        qc.ry(pi/4,0)             
        qc.ry(pi/2,1)  
        qc.ry(pi,2)             
        qc.cx (2,1)      
        qc.cx (1,0)
        qc.rz(pi/4,0)             
        qc.rz(pi/2,1)  
        qc.rz(pi,2)             
        qc.cx (2,1)      
        qc.cx (0,1)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==12:
        A = [[6, 5, 3], [5, 4, 5], [3, 5, 1]]
        qc = IQP(A)
        qc = qc.decompose ()
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==13:
        ansatz = EfficientSU2(3, entanglement='circular', reps=1)
        qc = QuantumCircuit(3)  # create a circuit and append the RY variational form
        qc.compose(ansatz, inplace=True)
        # parameterize
        param_dict = {}
        dividers = [-11,-6,15,10,7,20,-1,-20,-4,-10,17,-5]
        for key,divider in zip(qc.parameters,dividers):
            #for key in qc.parameters:
            #repeat = True
            #while repeat:
            #    divider = random.randint(-20,20)
            #    if divider != 0:
            #        repeat = False
            param_dict[key] = pi/divider
        qc = qc.bind_parameters(param_dict)
        qcd = qc.decompose()
        layers, num_layers = QCircuit_to_layers (qcd)
        return qcd, qcd.num_qubits, layers, num_layers
    elif ID==14:
        #qc = QuantumCircuit(3)
        qc = QuantumCircuit(4)
        # layer 0
        qc.h(0)
        qc.x(1)
        qc.h(2)
        qc.x(3)
        # layer 1
        qc.s(0)
        qc.s(1)
        qc.x(2)
        qc.h(3)
        # layer 2
        qc.h(0)
        qc.x(1)
        qc.s(2)
        qc.x(3)
        # layer 3
        qc.t(0)
        qc.t(1)
        qc.h(2)
        qc.t(3)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==15:
        theta = Parameter('θ')
        qc = QuantumCircuit(1)
        qc.rx(theta, 0)
        theta_val = pi / 4.
        qc = qc.bind_parameters({theta: theta_val})
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==16:
        qc = QuantumCircuit(3)
        qc.h(2)
        qc.cx(2,0)
        qc.cx(0,1)
        qc.cx(1,2)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==17:
        ansatz = EfficientSU2(4, entanglement='circular', reps=2)
        qc = QuantumCircuit(4)  # create a circuit and append the RY variational form
        qc.compose(ansatz, inplace=True)
        # parameterize
        param_dict = {}
        dividers = [10,-13,2,17,7,15,-8,-17,-7,-4,-4,1,-15,15,-4,14,5,9,12,-4,18,-6,11,18]
        for key,divider in zip(qc.parameters,dividers):
        #for key in qc.parameters:
        #    repeat = True
        #    while repeat:
        #        divider = random.randint(-20,20)
        #        if divider != 0:
        #            repeat = False
            param_dict[key] = pi/divider
        qc = qc.bind_parameters(param_dict)
        qcd = qc.decompose()
        layers, num_layers = QCircuit_to_layers (qcd)
        return qcd, qcd.num_qubits, layers, num_layers
    elif ID==18:
        A = [[1, 4, 6, 5, 3, 2], [4, 4, 5, 3, 7, 1], [6, 5, 2, 1, 4, 2], 
             [5, 3, 1, 3, 2, 4],[3, 7, 4, 2, 1, 5],[2, 1, 2, 4, 5, 7]]
        qc = IQP(A)
        qc = qc.decompose()
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==19:
        num_qubits = args[0]
        reps = args[1]
        ansatz = EfficientSU2(num_qubits, entanglement='linear', reps=reps)
        qc = QuantumCircuit(num_qubits)  # create a circuit and append the RY variational form
        qc.compose(ansatz, inplace=True)
        # parameterize
        param_dict = {}
        #dividers = [10,-13,2,17,7,15,-8,-17,-7,-4,-4,1,-15,15,-4,14,5,9,12,-4,18,-6,11,18]
        #for key,divider in zip(qc.parameters,dividers):
        random.seed (10000)
        for key in qc.parameters:
            repeat = True
            while repeat:
                divider = random.randint(-20,20)
                if divider != 0:
                    repeat = False
            param_dict[key] = pi/divider
        qc = qc.bind_parameters(param_dict)
        qcd = qc.decompose()
        layers, num_layers = QCircuit_to_layers (qcd)
        random.seed()
        return qcd, qcd.num_qubits, layers, num_layers
    elif ID==20:
        theta = Parameter('θ')
        qc = QuantumCircuit(2)
        qc.rx(theta, 0)
        qc.id(1)
        #qc.z(0)
        qc.ry(theta, 0)
        qc.z(1)
        qc.cx (0,1)
        theta_val = 7. *pi / 8.
        qc = qc.bind_parameters({theta: theta_val})
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==21:  #random circuit
        num_qubits = args[0]
        depth = args[1]
        qc = my_random_circuit(num_qubits=num_qubits, depth=depth, max_operands=2, measure=False, seed=393)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==23:  #random circuit with H, S, T and CNOT gates ; %of H is given as args[2]
        num_qubits = args[0]
        depth = args[1]
        Hpercent = args[2]
        qc = my_random_circuit_Cliff_T(num_qubits=num_qubits, depth=depth, max_operands=2, Hpercent= Hpercent, measure=False, seed=393)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==22:   # IQP
        # generate n x n interactions matrix. The powers of each T gate are given by the diagonal elements of the interactions matrix. The powers of the CS gates are given by the upper triangle of the interactions matrix.
        num_qubits = args[0]
        random.seed (10000)
        A = np.zeros((num_qubits,num_qubits), dtype=int)
        for i in range(num_qubits):
            for j in range (i, num_qubits):
                A[i][j] = random.randint(1,10)
                A[j][i] = A[i][j]
        qc = IQP(A)
        qc = qc.decompose ()
        layers, num_layers = QCircuit_to_layers (qc)
        random.seed()
        return qc, qc.num_qubits, layers, num_layers
    elif ID==100:   # test
        qc = circuit6(2)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==101:   # test 1 qubit: 3 layers Hadamard
        qc = QuantumCircuit(2)
        #Layer 0
        qc.h(0) 
        qc.h(1) 
        #Layer 1
        qc.h(0)
        qc.h(1) 
        #Layer 2
        qc.s(0) 
        qc.s(1)
        #Layer 3
        qc.h(0) 
        qc.h(1)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==102:   # Fake SU2
        qc = QuantumCircuit(2)
        #Layer 0
        qc.ry(pi/16,0)
        qc.ry(-pi,1)
        #Layer 1
        qc.rz(-pi/13,0)
        qc.rz(-pi/19,1)
        #Layer 2
        qc.cx(0,1)
        #Layer 3
        qc.ry(-pi/20,0)
        qc.ry(-pi/4,1)
        #Layer 4
        #qc.rz(-pi/9,0)
        #qc.rz(-pi/14,1)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID==103:   # Test P and CP
        qc = QuantumCircuit(2)
        #Layer 0
        qc.h(0)
        qc.cp(pi/3, 0, 1)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    
    elif ID==300:   # IQP inversion test (for kernel methods: computes the overlap as the probability of)
        # generate n x n interactions matrix. The powers of each T gate are given by the diagonal elements of the interactions matrix. The powers of the CS gates are given by the upper triangle of the interactions matrix.
        num_qubits = args[0]
        random.seed (10000)
        A = np.zeros((num_qubits,num_qubits), dtype=int)
        for i in range(num_qubits):
            for j in range (i, num_qubits):
                A[i][j] = random.randint(1,10)
                A[j][i] = A[i][j]
        qcA = IQP(A)
        B = np.zeros((num_qubits,num_qubits), dtype=int)
        for i in range(num_qubits):
            for j in range (i, num_qubits):
                B[i][j] = random.randint(1,10)
                B[j][i] = B[i][j]
        qcB = IQP(B)
        qcB_inv = qcA.inverse()
        qc = qcB.compose(qcB_inv)
        qc = qc.decompose ()
        qc = _remove_H_pairs (qc)
        layers, num_layers = QCircuit_to_layers (qc)
        random.seed()
        return qc, qc.num_qubits, layers, num_layers
    elif ID == 400:   # HSA
        n = args[0]  # number of qubits
        # n must be even
        n -= 0 if n%2==0 else 1
        n_ccz = args[1]    # implies t_count = 14 * n_ccz
        g = args[2] # number of Z and CZ gates in the {Z,CZ}-layers in-between Toffoli gates 
        # Generating a random hidden string
        s = []
        for i in range(0,n,2):
            s.append(0)
            s.append(1)
        print (s)
        random.seed (10000)            
        qc = HSA_gen (n, n_ccz, g, s, Toffoli= False, measure=False, remove_pairs_H=True)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID == 4100:   # HSA
        n = args[0]  # number of qubits
        # n must be even
        n -= 0 if n%2==0 else 1
        n_ccz = args[1]    # implies t_count = 14 * n_ccz
        g = args[2] # number of Z and CZ gates in the {Z,CZ}-layers in-between Toffoli gates 
        Y = args[3]
        # Generating a random hidden string
        s = [0] * n
        i = 0
        while Y>0:
            s[i] = Y%2
            Y = Y //2
            i= i+1
        print (s)
        random.seed (10000)            
        qc = HSA_gen (n, n_ccz, g, s, Toffoli= False, measure=False, remove_pairs_H=True)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
    elif ID == 500:
        n = args[0]  # number of qubits
        nCycles = args[1]    # layers of CZ
        random.seed (10000)            
        qc = ran_gen (n, nCycles, measure=False, remove_pairs_H=True)
        layers, num_layers = QCircuit_to_layers (qc)
        return qc, qc.num_qubits, layers, num_layers
        
###################
#
# My random circuit generator
#
# Adapted from qiskit.circuit.random.random_circuit
#

import numpy as np

from qiskit.circuit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit.circuit.library.standard_gates import (
    IGate,
    XGate,
    ZGate,
    HGate,
    SGate,
    TGate,
    RXGate,
    RYGate,
    RZGate,
    CXGate,
    CZGate
)
from qiskit.circuit.exceptions import CircuitError


def my_random_circuit(
    num_qubits, depth, max_operands=2, measure=False, seed=None
):
    """Generate random circuit of arbitrary size and form.

    This function will generate a random circuit by randomly selecting gates
    from the set of standard gates in :mod:`qiskit.extensions`. For example:

    .. jupyter-execute::

        from qiskit.circuit.random import random_circuit

        circ = random_circuit(2, 2, measure=True)
        circ.draw(output='mpl')

    Args:
        num_qubits (int): number of quantum wires
        depth (int): layers of operations (i.e. critical path length)
        max_operands (int): maximum operands of each gate (between 1 and 2)
        measure (bool): if True, measure all qubits at the end
        seed (int): sets random seed (optional)

    Returns:
        QuantumCircuit: constructed circuit

    Raises:
        CircuitError: when invalid options given
    """
    if max_operands < 1 or max_operands > 2:
        raise CircuitError("max_operands must be between 1 and 2")

    one_q_ops = [
        IGate,
        XGate,
        ZGate,
        HGate,
        SGate,
        TGate,
        RXGate,
        RYGate,
        RZGate,
    ]
    one_param = [RXGate, RYGate, RZGate]
    two_q_ops = [CXGate, CZGate]

    qr = QuantumRegister(num_qubits, "q")
    qc = QuantumCircuit(num_qubits)

    if measure:
        cr = ClassicalRegister(num_qubits, "c")
        qc.add_register(cr)

    if seed is None:
        seed = np.random.randint(0, np.iinfo(np.int32).max)
    rng = np.random.default_rng(seed)

    # apply arbitrary random operations at every depth
    for _ in range(depth):
        # choose either 1 or 2 qubits for the operation
        remaining_qubits = list(range(num_qubits))
        rng.shuffle(remaining_qubits)
        while remaining_qubits:
            max_possible_operands = min(len(remaining_qubits), max_operands)
            possible_operands = range(max_possible_operands)
            len_possible_operands = len (possible_operands)
            num_operands = rng.choice(possible_operands, p=[0.8,0.2] if len_possible_operands==2 else [1.]) + 1
            operands = [remaining_qubits.pop() for _ in range(num_operands)]
            if num_operands == 1:
                operation = rng.choice(one_q_ops)
            elif num_operands == 2:
                operation = rng.choice(two_q_ops)
            if operation in one_param:
                num_angles = 1
            else:
                num_angles = 0
            angles = [rng.uniform(0, 2 * np.pi) for x in range(num_angles)]
            register_operands = [qr[i] for i in operands]
            op = operation(*angles)

            qc.append(op, register_operands)

    if measure:
        qc.measure(qr, cr)

    return qc

def my_random_circuit_Cliff_T(
    num_qubits, depth, max_operands=2, Hpercent=0.3, measure=False, seed=None
):
    """Generate random circuit of arbitrary size and form.

    This function will generate a random circuit by randomly selecting gates
    from the set of standard gates in :mod:`qiskit.extensions`. For example:

    .. jupyter-execute::

        from qiskit.circuit.random import random_circuit

        circ = random_circuit(2, 2, measure=True)
        circ.draw(output='mpl')

    Args:
        num_qubits (int): number of quantum wires
        depth (int): layers of operations (i.e. critical path length)
        max_operands (int): maximum operands of each gate (between 1 and 2)
        Hpercent (float): percentage of Hadamards (branching)
        measure (bool): if True, measure all qubits at the end
        seed (int): sets random seed (optional)

    Returns:
        QuantumCircuit: constructed circuit

    Raises:
        CircuitError: when invalid options given
    """
    if max_operands < 1 or max_operands > 2:
        raise CircuitError("max_operands must be between 1 and 2")

    ops = [
        HGate,
        SGate,
        TGate,
        CXGate
    ]
    
    ops_ndx = [0, 1, 2, 3]
    
    notHpercent = (1.-Hpercent) / 3.
    select_ops_prob = [Hpercent, notHpercent, notHpercent, notHpercent]

    notHpercent_not2ops = (1.-Hpercent) / 2.
    select_ops_prob_not2ops = [Hpercent, notHpercent_not2ops, notHpercent_not2ops, 0.]

    qr = QuantumRegister(num_qubits, "q")
    qc = QuantumCircuit(num_qubits)

    if measure:
        cr = ClassicalRegister(num_qubits, "c")
        qc.add_register(cr)

    if seed is None:
        seed = np.random.randint(0, np.iinfo(np.int32).max)
    rng = np.random.default_rng(seed)

    # apply arbitrary random operations at every depth
    for _ in range(depth):
        # choose either 1 or 2 qubits for the operation
        remaining_qubits = list(range(num_qubits))
        rng.shuffle(remaining_qubits)
        while remaining_qubits:
            # select the gate
            max_possible_operands = min(len(remaining_qubits), max_operands)
            operation_ndx = rng.choice(ops_ndx, p=select_ops_prob if max_possible_operands > 1 else select_ops_prob_not2ops)
            operation = ops[operation_ndx]
            num_operands = 1 if operation_ndx <3 else 2
            operands = [remaining_qubits.pop() for _ in range(num_operands)]
            register_operands = [qr[i] for i in operands]

            qc.append(operation(), register_operands)

    if measure:
        qc.measure(qr, cr)

    return qc
    