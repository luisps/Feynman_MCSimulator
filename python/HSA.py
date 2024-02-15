import qiskit
import random

def oracle_g(n, n_ccz, g, Toffoli):
    """Oragel O_g

    Args:
        n (int): Number of qubits in the quantum circuit
        n_ccz (int): Number of Toffoli (or CZZ) gates
        g (int): Number of Z and CZ gates in the {Z,CZ}-layers

    Returns:
        QuantumCircuit: The O_g oracle quantum circuit
    """
    qr_g = qiskit.QuantumRegister(int(n / 2))

    og = qiskit.QuantumCircuit(qr_g, name='oracle_g')

    gates = ['z', 'cz']
    qubits = [i for i in range(int(n / 2))]
    for i in range(g):
        gate = random.sample(gates, 1)
        if gate[0] == 'z':
            target = random.sample(qubits, 1)
            og.s(qr_g[target])
            og.s(qr_g[target])
        elif gate[0] == 'cz':
            control, target = random.sample(qubits, 2)
            og.h(qr_g[target])
            og.cx(qr_g[control], qr_g[target])
            og.h(qr_g[target])

    for i in range(n_ccz):
        c1, c2, t = random.sample(qubits, 3)
        if Toffoli:
            og.h(qr_g[t])
            og.ccx(qr_g[c1], qr_g[c2], qr_g[t])
            og.h(qr_g[t])
        else:
            # cx from c2 to target
            og.cx (qr_g[c2], qr_g[t])
            # t dagger on target
            og.t (qr_g[t])
            og.s (qr_g[t])
            og.s (qr_g[t])
            og.s (qr_g[t])
            # cx from c1 to target
            og.cx (qr_g[c1], qr_g[t])
            # t on target
            og.t (qr_g[t])
            # cx from c2 to target
            og.cx (qr_g[c2], qr_g[t])
            # t dagger on target
            og.t (qr_g[t])
            og.s (qr_g[t])
            og.s (qr_g[t])
            og.s (qr_g[t])
            # cx from c1 to target
            og.cx (qr_g[c1], qr_g[t])
            # t on target and c2
            og.t (qr_g[t])
            og.t (qr_g[c2])
            # cx from c1 to c2
            og.cx (qr_g[c1], qr_g[c2])
            # t dagger on c2
            og.t (qr_g[c2])
            og.s (qr_g[c2])
            og.s (qr_g[c2])
            og.s (qr_g[c2])
            # t on c1
            og.t (qr_g[c1])
            # cx from c1 to c2
            og.cx (qr_g[c1], qr_g[c2])            
        for i in range(g):
            gate = random.sample(gates, 1)
            if gate[0] == 'z':
                target = random.sample(qubits, 1)
                og.s(qr_g[target])
                og.s(qr_g[target])
            elif gate[0] == 'cz':
                control, target = random.sample(qubits, 2)
                og.h(qr_g[target])
                og.cx(qr_g[control], qr_g[target])
                og.h(qr_g[target])

    return og


def oracle_f(n, og):
    """Oracle O_f

    Args:
        n (int): The number of qubits in the quantum circuit
        og (QuantumCircuit): The O_g oracle quantum circuit

    Returns:
        QuantumCircuit: The O_f oracle quantum circuit
    """
    qr_f = qiskit.QuantumRegister(n)

    of = qiskit.QuantumCircuit(qr_f, name='oracle_f')

    of = of.compose(og, qubits=qr_f[0:int(n / 2)])  # compose

    for i in range(int(n / 2)):
        of.h(qr_f[i + int(n / 2)])
        of.cx(qr_f[i], qr_f[i + int(n / 2)])
        of.h(qr_f[i + int(n / 2)])

    return of


def oracle_fp(n, og, s):
    """Oracle O_f^{\prime}

    Args:
        n (int): The number of qubits in the quantum circuit
        og (QuantumCircuit): The O_g oracle quantum circuit

    Returns:
        QuantumCircuit: The O_f^{\prime} oracle quantum circuit
    """
    qr_fp = qiskit.QuantumRegister(n)

    ofp = qiskit.QuantumCircuit(qr_fp, name='oracle_fp')

    for i in range(n):
        if s[i] == 1:
            ofp.s(qr_fp[i])
            ofp.s(qr_fp[i])

    ofp = ofp.compose(og, qubits=qr_fp[int(n / 2):n])  # compose

    for i in range(int(n / 2)):
        ofp.h(qr_fp[i + int(n / 2)])
        ofp.cx(qr_fp[i], qr_fp[i + int(n / 2)])
        ofp.h(qr_fp[i + int(n / 2)])

    return ofp

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

def HSA_gen (n, n_ccz, g, s, Toffoli = True,measure = True, remove_pairs_H=False):
    og = oracle_g(n, n_ccz, g, Toffoli)
    of = oracle_f(n, og)
    ofp = oracle_fp(n, og, s)
        
    qr = qiskit.QuantumRegister(n, 'qr')
    if measure:
        cr = qiskit.ClassicalRegister(n, 'cr')
        qc = qiskit.QuantumCircuit(qr, cr, name='U')
    else:
        qc = qiskit.QuantumCircuit(qr, name='U')
        
    for i in range(n):
        qc.h(i)
        
    qc = qc.compose(of, qr[0: n])
        
    for i in range(n):
        qc.h(i)
        
    qc = qc.compose(ofp, qr[0: n])
        
    for i in range(n):
        qc.h(i)
        if measure:
            qc.measure(i, n - 1 - i)

    if remove_pairs_H:
        qc = remove_H_pairs (qc)
        
    return qc