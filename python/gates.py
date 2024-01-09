
# methods to sample I, H, S, Z, X, T, CZ, CX, CP, RX, RY, RZ, P gates


# For each type of gate the <gate>_w_prob () method (e.g H_w_prob () )
# takes the following parameters:
#     . i_qubits: TUPLE with the state of the input qubits, i.e, 0 or 1
#     . o_qubits: TUPLE with the state of the output qubits, i.e, 0 or 1
#     . distribution: identifier of the distribution to use to sample the gate
#                     currently only 'uniform' is supported
#
# The <gate>_w_prob() method returns two values:
#     . w - the amplitude of the transition from i_qubits to o_qubits for that gate
#     . prob - the probability of the above transition given the used distribution

# Definitions of the transition weights matrix for each gate

from math import sqrt, pi, cos, sin
import cmath
import numpy as np

####################################
# Unitaries for the supported gates
#

# Hadamard
REC_SQRT_2 = 1. / sqrt(2) 
H = [[REC_SQRT_2, REC_SQRT_2],[REC_SQRT_2, -REC_SQRT_2]]

# S 
S = [[1., 0.], [0., complex(0,1)]]

# T 
T = [[1., 0.], [0., cmath.exp(complex(0,pi/4.))]]

# I
I = [[1., 0.], [0., 1.]]

# Z 
Z = [[1., 0.], [0., -1.]]

# X 
X = [[0., 1.], [1., 0.]]

# CX (control qubit,target qubit)
CX = [[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 0., 1.], [0., 0., 1., 0.]]

# CZ 
CZ = [[1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., -1.]]



# gate_evaluate_w (gate_name, input_qbs, output_qbs, , params=[], unitary=[]):
# given:
#     . the gate name
#     . a tuple with the input qubits
#     . a tuple with the output qubits
#
# return the amplitude w associated with transitioning from input to output
def gate_evaluate_w (gate_name, input_qbs, output_qbs, params=[], unitary=[]):
    if gate_name =='h':
        w,_ = H_w_prob(input_qbs, output_qbs)
    elif gate_name =='id':
        w,_ = I_w_prob(input_qbs, output_qbs)
    elif gate_name =='s':
        w,_ = S_w_prob(input_qbs, output_qbs)
    elif gate_name =='z':
        w,_ = Z_w_prob(input_qbs, output_qbs)
    elif gate_name =='x':
        w,_ = X_w_prob(input_qbs, output_qbs)
    elif gate_name =='t':
        w,_ = T_w_prob(input_qbs, output_qbs)
    elif gate_name =='cx':
        w,_ = CX_w_prob(input_qbs, output_qbs)
    elif gate_name =='cz':
        w,_ = CZ_w_prob(input_qbs, output_qbs)
    elif gate_name =='cp':
        w,_ = CP_w_prob(input_qbs, output_qbs, unitary=unitary)       
    elif gate_name =='rx':
        w,_ = RX_w_prob(input_qbs, output_qbs, unitary=unitary)
    elif gate_name =='ry':
        w,_ = RY_w_prob(input_qbs, output_qbs, unitary=unitary)
    elif gate_name =='rz':
        w,_ = RZ_w_prob(input_qbs, output_qbs, unitary=unitary)       
    elif gate_name =='p':
        w,_ = P_w_prob(input_qbs, output_qbs, unitary=unitary)       
    else:  # non defined gate: raise error
        raise ValueError('Undefined gate name = {0}'.format(gate_name))
    return w

# gate_evaluate_w_prob (gate_name, input_qbs, output_qbs, , params=[], unitary=[]):
# given:
#     . the gate name
#     . a tuple with the input qubits
#     . a tuple with the output qubits
#
# return the amplitude w and the prob associated with transitioning from input to output
def gate_evaluate_w_prob (gate_name, input_qbs, output_qbs, params=[], unitary=[], pdfs=None):
    if gate_name =='h':
        w,prob = H_w_prob(input_qbs, output_qbs)
    elif gate_name =='id':
        w,prob = I_w_prob(input_qbs, output_qbs)
    elif gate_name =='s':
        w,prob = S_w_prob(input_qbs, output_qbs)
    elif gate_name =='z':
        w,prob = Z_w_prob(input_qbs, output_qbs)
    elif gate_name =='x':
        w,prob = X_w_prob(input_qbs, output_qbs)
    elif gate_name =='t':
        w,prob = T_w_prob(input_qbs, output_qbs)
    elif gate_name =='cx':
        w,prob = CX_w_prob(input_qbs, output_qbs)
    elif gate_name =='cz':
        w,prob = CZ_w_prob(input_qbs, output_qbs)
    elif gate_name =='cp':
        w,prob = CP_w_prob(input_qbs, output_qbs, unitary=unitary)       
    elif gate_name =='rx':
        w,prob = RX_w_prob(input_qbs, output_qbs, unitary=unitary, pdfs=pdfs)
    elif gate_name =='ry':
        w,prob = RY_w_prob(input_qbs, output_qbs, unitary=unitary, pdfs=pdfs)
    elif gate_name =='rz':
        w,prob = RZ_w_prob(input_qbs, output_qbs, unitary=unitary, pdfs=pdfs)       
    elif gate_name =='p':
        w,prob = P_w_prob(input_qbs, output_qbs, unitary=unitary)       
    else:  # non defined gate: raise error
        raise ValueError('Undefined gate name = {0}'.format(gate_name))
    return w, prob

# gate_evaluate_w_sample (gate_name, input_qbs, rnd=None, unitary=[], pdfs=[], smpl_dir="DIRECT"):
# given:
#     . the gate name
#     . a tuple with the input qubits
#     . a random number in [0,1[
#     . the gate unitary
#     . the gate pdfs
#.    . the direction of sampling (DIRECT, REVERSE)
#
# return:
#     . a tuple output_qbs with the output qubits
#     . the amplitude w associated with transtioning from input to output
#     . the probability prob of this transition happening
def gate_evaluate_w_sample (gate_name, input_qbs, rnd=None, unitary=[], pdfs=[], smpl_dir="DIRECT"):
    if gate_name =='h':          # invariant with the sampling direction
        output_qbs, w, prob = H_w_sample(input_qbs, rnd=rnd)
    elif gate_name =='id':
        output_qbs = input_qbs    # invariant with the sampling direction
        w = prob = 1.
    elif gate_name =='s':         # invariant with the sampling direction
        output_qbs, w, prob = S_w_sample(input_qbs)
    elif gate_name =='z':         # invariant with the sampling direction
        output_qbs, w, prob = Z_w_sample(input_qbs)
    elif gate_name =='x':         # invariant with the sampling direction
        output_qbs, w, prob = X_w_sample(input_qbs)
    elif gate_name =='t':         # invariant with the sampling direction
        output_qbs, w, prob = T_w_sample(input_qbs)
    elif gate_name =='cx':        # invariant with the sampling direction
        output_qbs, w, prob = CX_w_sample(input_qbs)
    elif gate_name =='cz':        # invariant with the sampling direction
        output_qbs, w, prob = CZ_w_sample(input_qbs)
    elif gate_name =='cp':        # invariant with the sampling direction
        output_qbs, w, prob = CP_w_sample(input_qbs, unitary=unitary)       
    elif gate_name in ['rx','ry','rz']:
        output_qbs, w, prob = rotation_w_sample(input_qbs, rnd=rnd, unitary=unitary, pdfs=pdfs, smpl_dir=smpl_dir)
    elif gate_name =='p':         # invariant with the sampling direction
        output_qbs, w, prob = P_w_sample(input_qbs, unitary=unitary)       
    else:  # non defined gate: raise error
        raise ValueError('Undefined gate name = {0}'.format(gate_name))
    return output_qbs, w, prob


def X_w_sample(input_qbs):
    output_qbs = (1-input_qbs[0],)
    w = prob = 1.
    return output_qbs, w, prob

def Z_w_sample(input_qbs):
    output_qbs = input_qbs
    prob = 1.   # no branching gate
    w = Z[input_qbs[0]][input_qbs[0]]  # diagonal gate
    return output_qbs, w, prob

def S_w_sample(input_qbs):
    output_qbs = input_qbs
    prob = 1.   # no branching gate
    w = S[input_qbs[0]][input_qbs[0]]  # diagonal gate
    return output_qbs, w, prob

def T_w_sample(input_qbs):
    output_qbs = input_qbs
    prob = 1.   # no branching gate
    w = T[input_qbs[0]][input_qbs[0]]  # diagonal gate
    return output_qbs, w, prob

def H_w_sample (input_qbs, rnd=None):
    # Stochastically select the output qubit value
    output_qbs = (1,) if rnd>= 0.5 else (0,)
    # the output qubit selects the row
    row = output_qbs[0]
    # the input qubit selects the column
    col = input_qbs[0]    
    # the weight of the transition is given by the matrix
    w = H[row][col]
    # the probability for the Hadamard gate and 'uniform' is 0.5
    prob = 0.5
    return output_qbs, w, prob

def CX_w_sample (input_qbs):
    # input_qbs[0] is the control qubit
    # input_qbs[1] is the target qubit
    o = input_qbs[1] if input_qbs[0]==0 else 1-input_qbs[1]
    output_qbs = (input_qbs[0], o)
    w = 1.
    prob = 1.
    return output_qbs, w, prob
# end CX_w_sample

def CZ_w_sample (input_qbs):
    # input_qbs[0] is the control qubit
    # input_qbs[1] is the target qubit
    index = input_qbs[0]*2 + input_qbs[1]
    output_qbs = input_qbs  # no branching gate
    w = CZ[index][index]    # diagonal gate
    prob = 1.
    return output_qbs, w, prob
# end CZ_w_sample

def P_w_sample(input_qbs, unitary=[]):
    output_qbs = input_qbs
    prob = 1.   # no branching gate
    w = unitary[input_qbs[0]][input_qbs[0]]  # diagonal gate
    return output_qbs, w, prob

def CP_w_sample (input_qbs, unitary=[]):
    # the input qubits (2 of them) select the column
    # input_qbs[0] is the control qubit: most significant in the 2 bits string
    # input_qbs[1] is the target qubit: least significant in the 2 bits string
    index = input_qbs[0]*2 + input_qbs[1]
    output_qbs = input_qbs      # no branching gate
    # the weight of the transition is given by the matrix
    w = unitary[index][index]   # diagonal gate
    prob = 1.0
    return output_qbs, w, prob
# end CP_w_sample

def rotation_w_sample (input_qbs, rnd=None, unitary=[], pdfs=[], smpl_dir="DIRECT"):
    # Remember that the input qubit selects which column of pdfs to use
    col = input_qbs[0]    
    # Stochastically generate the output qubit
    # the output qubit selects the row
    output_qbs = (1,) if rnd>= pdfs[0][col] else (0,)
    # the output qubit selects the row
    row = output_qbs[0]
    # RX matrix:
    # [ cos(half_theta)     -i sin(half_theta)  ]
    # [ -i sin(half_theta)  cos(half_theta)     ]
    # the weight of the transition is given by the matrix

    # sampling direction
    # Probabilities do not change with the sampling direction
    # Amplitudes do not change
    # but this method parameterization requires exchanging input with output if "REVERSE"
    w = unitary[row][col] if smpl_dir == "DIRECT" else unitary[col][row]
    prob = pdfs[row][col]
    
    return output_qbs, w, prob
# end rotation_w_sample

def H_w_prob (i_qubits, o_qubits):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # the weight of the transition is given by the matrix
    w = H[row][col]
    # the probability for the Hadamard gate and 'uniform' is 0.5
    prob = 0.5
    return w, prob
# end H_w_prob

def S_w_prob (i_qubits, o_qubits):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # the weight of the transition is given by the matrix
    w = S[row][col]
    # the probability for the S gate and 'uniform' is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.0
    return w, prob
# end S_w_prob

def Z_w_prob (i_qubits, o_qubits):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # the weight of the transition is given by the matrix
    w = Z[row][col]
    # the probability for the Z gate and 'uniform' is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.0
    return w, prob
# end Z_w_prob

def X_w_prob (i_qubits, o_qubits):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # the weight of the transition is given by the matrix
    w = X[row][col]
    # the probability for the X gate and 'uniform' is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.0
    return w, prob
# end X_w_prob

def P_w_prob (i_qubits, o_qubits, params= [], unitary=[]):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # P matrix:
    # [[1., 0.] , 
    # [0., complex(cos(psi), sin(psi))]]
    # the weight of the transition is given by the matrix
    w = unitary[row][col]
    # the probability for the P gate and 'uniform' is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.
    return w, prob
# end P_w_prob

def RX_w_prob (i_qubits, o_qubits, params= [], unitary=[], pdfs=None):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # RX matrix:
    # [ cos(half_theta)     -i sin(half_theta)  ]
    # [ -i sin(half_theta)  cos(half_theta)     ]
    # the weight of the transition is given by the matrix
    w = unitary[row][col]
    # the probability for the RX gate is :
    prob = pdfs[row][col] if not pdfs is None else 1.
    return w, prob
# end RX_w_prob

def RY_w_prob (i_qubits, o_qubits, params= [], unitary=[], pdfs=None):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # RY matrix:
    # [[cos(half_theta), -sin(half_theta)] , 
    #  [sin(half_theta), cos(half_theta)]]
    # the weight of the transition is given by the matrix
    w = unitary[row][col]
    # the probability for the RY gate is :
    prob = pdfs[row][col] if not pdfs is None else 1.
    return w, prob
# end RY_w_prob

def RZ_w_prob (i_qubits, o_qubits, params= [], unitary=[], pdfs=None):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # RZ matrix:
    # [[complex(cos(half_theta), -sin(half_theta)), 0.] , 
    #  [0., complex(cos(half_theta),  sin(half_theta))]]
    # the weight of the transition is given by the matrix
    w = unitary[row][col]
    # the probability for the RZ gate is :
    prob = pdfs[row][col] if not pdfs is None else 1.
    return w, prob
# end RZ_w_prob

def T_w_prob (i_qubits, o_qubits):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # the weight of the transition is given by the matrix
    w = T[row][col]
    # the probability for the T gate and 'uniform' is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.0
    return w, prob
# end T_w_prob

def I_w_prob (i_qubits, o_qubits):
    # the output qubit selects the row
    row = o_qubits[0]
    # the input qubit selects the column
    col = i_qubits[0]    
    # the weight of the transition is given by the matrix
    w = I[row][col]
    # the probability for the I gate and 'uniform' is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.0
    return w, prob
# end I_w_prob

def CP_w_prob (i_qubits, o_qubits, params= [], unitary=[]):
    # the output qubits select the row
    # o_qubits[0] is the control qubit: most significant in the 2 bits string
    # o_qubits[1] is the target qubit: least significant in the 2 bits string
    row = o_qubits[0]*2 + o_qubits[1]
    # the input qubits (2 of them) select the column
    # i_qubits[0] is the control qubit: most significant in the 2 bits string
    # i_qubits[1] is the target qubit: least significant in the 2 bits string
    col = i_qubits[0]*2 + i_qubits[1]    
    # the weight of the transition is given by the matrix
    w = unitary[row][col]
    # the probability for the CP gate is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.0
    return w, prob
# end CP_w_prob

def CX_w_prob (i_qubits, o_qubits):
    # the output qubits select the row
    # o_qubits[0] is the control qubit: most significant in the 2 bits string
    # o_qubits[1] is the target qubit: least significant in the 2 bits string
    row = o_qubits[0]*2 + o_qubits[1]
    # the input qubits (2 of them) select the column
    # i_qubits[0] is the control qubit: most significant in the 2 bits string
    # i_qubits[1] is the target qubit: least significant in the 2 bits string
    col = i_qubits[0]*2 + i_qubits[1]    
    # the weight of the transition is given by the matrix
    w = CX[row][col]
    # the probability for the CX gate is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.0
    return w, prob
# end CX_w_prob

def CZ_w_prob (i_qubits, o_qubits):
    # the output qubits (2 of them) select the row
    row = o_qubits[1]*2 + o_qubits[0]
    # the input qubits (2 of them) select the column
    col = i_qubits[1]*2 + i_qubits[0]    
    # the weight of the transition is given by the matrix
    w = CZ[row][col]
    # the probability for the CZ gate is :
    #   if w = 0. then prob = 0.0
    #   if w != 0. then prob = 1.0
    prob = 0.0 if w==0.0 else 1.0
    return w, prob
# end CZ_w_prob

