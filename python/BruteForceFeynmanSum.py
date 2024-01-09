from hash_paths import my_hash_inverse, my_hash 
from gates import gate_evaluate_w

import csv

def BruteForce (num_qubits, num_layers, layers, initial_bin_str, final_bin_str, P, N):
    __DEBUG__ = False
    format_str = '{{0:0{0}b}}'.format(num_qubits)

    amplitude = complex (0., 0.)
    non_zero_paths_ndx = []

    # create a list of lists to hold the values of the qubits
    # at each intermediate state
    qb_v = []
    for qb in range(num_qubits):
        qb_v.append([])
        for l in range(num_layers+1):
            qb_v[qb].append(None)

    for qb in range(num_qubits):
        #inital state
        qb_v[num_qubits-qb-1][0] = int(initial_bin_str[qb])
        #final state
        qb_v[num_qubits-qb-1][num_layers] = int(final_bin_str[qb])

    # for depicting how amplitudes are distributed across paths
    amplitudes = []
    # iterate over all paths
    for path_ndx in range(P):
    #for path_ndx in range(1):
        path = my_hash_inverse (path_ndx, N, num_layers-1)
        if __DEBUG__:
            print ('{0}: {1}'.format(path_ndx, path))
            print (initial_bin_str)
        for l in range(1,num_layers):
            bstring = format_str.format(path[l-1])
            if __DEBUG__:
                print (bstring)
            for qb in range(num_qubits):
                qb_v[qb][l] = int(bstring[num_qubits-qb-1])
        if __DEBUG__:
            print (final_bin_str)
        
        if __DEBUG__:   # print the states bit strings as stored in qb_v
            for l in range(0,num_layers+1):
                print ('state {0}'.format(l))
                for qb in range(num_qubits):
                    print (qb_v[qb][l], end = ' ')
                print ()

        # compute the product across layers for the current path
        path_throughput = complex(1., 0.)  # path initially has weight 1. + 0j
        for l in range(num_layers):
            if __DEBUG__:
                print ('Layer {0}'.format(l))
            # compute the product across gates for the current path and layer
            layer_weight = complex(1., 0.)  # layer initially has weight 1
        
            layer = layers[l]
        
            processed_qubits = []   # some gates have 2 qubits: 
                                    # make sure we don't process twice
            # iterate over all qubits 
            for qubit in range(num_qubits):
                                
                if qubit in processed_qubits:
                    continue
                
                gate = layer[qubit]
                gate_name = gate[0]
                qubits = gate[1]  # list of addressed qubits
                params = gate[2]  # list of parameters to the gate
                unitary = gate[3] # unitary for parameterized gates
            
                if __DEBUG__:
                    print ('Qubit {0}, gate={1}'.format(qubit, gate), end='; ')

                for qb in qubits:
                    processed_qubits.append(qb)
                        
                # 2 qubit gates
                if gate_name in ['cx','cz','cp']:
                    input_qbs = (qb_v[qubits[0]][l], qb_v[qubits[1]][l])
                    output_qbs = (qb_v[qubits[0]][l+1],qb_v[qubits[1]][l+1])
                else: 
                    input_qbs = (qb_v[qubit][l], )
                    output_qbs = (qb_v[qubit][l+1],)
                w = gate_evaluate_w (gate_name, input_qbs, output_qbs, 
                                     params=params, unitary=unitary)

                if __DEBUG__:
                    print ('Gate {0} with input={1} and output={2} has amplitude={3}'.format(gate_name, input_qbs, output_qbs, w))
                
                if w==0.:
                    layer_weight = 0.
                    path_throughput = 0.
                    break
                else:
                    layer_weight *= w
                                    
            # end iterating over the layer qubits
            if __DEBUG__:
                print ()
                
            if layer_weight == 0.:
                path_throughput = 0.
                break
            else:
                path_throughput *= layer_weight
            
        # end iterating over layers
        #print ('{0:.5}'.format(path_throughput))
        amplitudes.append(path_throughput)
        if path_throughput != 0.:
            amplitude += path_throughput
            non_zero_paths_ndx.append(path_ndx)

            if __DEBUG__:
                print (path_ndx, path)
                print (path_throughput)
                print ()

    # end iterating over paths
    
    return amplitude, amplitudes, non_zero_paths_ndx

def BruteForcePDF (num_qubits, num_layers, layers, initial_bin_str, final_bin_str, P, N):
    __DEBUG__ = False
    format_str = '{{0:0{0}b}}'.format(num_qubits)

    amplitude = complex (0., 0.)
    non_zero_paths_ndx = []
    EPSILON = 1e-5

    # create a LIST OF LISTS to hold the values of the qubits
    # at each intermediate state
    qb_v = []
    for qb in range(num_qubits):
        qb_v.append([])
        for l in range(num_layers+1):
            qb_v[qb].append(None)

    # fill in initial and final state (qb_v[][0] and qb_v[][num_layers]
    for qb in range(num_qubits):
        #inital state
        qb_v[num_qubits-qb-1][0] = int(initial_bin_str[qb])
        #final state
        qb_v[num_qubits-qb-1][num_layers] = int(final_bin_str[qb])

    # for storing amplitudes for all paths
    amplitudes = []
    # for storing pdfs for all paths
    pdfs = []
    # iterate over all paths
    for path_ndx in range(P):
    #for path_ndx in range(1):
        if (path_ndx & 0x03FF)==0: print ('.', end=' ')
        path = my_hash_inverse (path_ndx, N, num_layers-1)
        if __DEBUG__:
            print ('{0}: {1}'.format(path_ndx, path))
            print (initial_bin_str)
        # fill in qb_v[0..num_qubits-1] with 0 or 1 for all intermediate states
        for l in range(1,num_layers):
            bstring = format_str.format(path[l-1])
            if __DEBUG__:
                print (bstring)
            for qb in range(num_qubits):
                qb_v[qb][l] = int(bstring[num_qubits-qb-1])
        if __DEBUG__:
            print (final_bin_str)
        
        if __DEBUG__:   # print the states bit strings as stored in qb_v
            for l in range(0,num_layers+1):
                print ('state {0}'.format(l))
                for qb in range(num_qubits):
                    print (qb_v[qb][l], end = ' ')
                print ()

        # compute the product across layers for the current path
        path_throughput = complex(1., 0.)  # path initially has weight 1+0j
        path_pdf = 1.   # path initially has pdf 1.
        for l in range(num_layers):
            if __DEBUG__:
                print ('Layer {0}'.format(l))
            # compute the product across gates for the current path and layer
            layer_weight = complex(1., 0.)  # layer initially has weight 1
            layer_pdf = 1.  # layer initially has pdf 1.
        
            layer = layers[l]
        
            processed_qubits = []   # some gates have 2 qubits: 
                                    # make sure we don't process twice
            # iterate over all qubits 
            for qubit in range(num_qubits):
                                
                if qubit in processed_qubits:
                    continue
                
                gate = layer[qubit]
                gate_name = gate[0]
                qubits = gate[1]  # list of addressed qubits
                params = gate[2]  # list of parameters to the gate
                unitary = gate[3] # unitary for parameterized gates
            
                if __DEBUG__:
                    print ('Qubit {0}, gate={1}'.format(qubit, gate), end='; ')

                for qb in qubits:
                    processed_qubits.append(qb)
                        
                # 2 qubit gates
                if gate_name in ['cx','cz','cp']:
                    input_qbs = (qb_v[qubits[0]][l], qb_v[qubits[1]][l])
                    output_qbs = (qb_v[qubits[0]][l+1],qb_v[qubits[1]][l+1])
                else: 
                    input_qbs = (qb_v[qubit][l], )
                    output_qbs = (qb_v[qubit][l+1],)
                w = gate_evaluate_w (gate_name, input_qbs, output_qbs, 
                                     params=params, unitary=unitary)
                
                # the probability of selecting this transition is:
                #   ... the square of the absolute value of the amplitude if this is not the last state
                #   ... 1 if we are addressing the last state, since this is fixed
                pdf = 1. if l==num_layers-1 else abs(w)**2

                if __DEBUG__:
                    print ('Gate {0} with input={1} and output={2} has amplitude={3}'.format(gate_name, input_qbs, output_qbs, w))
                
                layer_weight *= w
                layer_pdf *= pdf
                
                if pdf==0.: break
                                    
            # end iterating over the layer qubits

            path_throughput *= layer_weight
            path_pdf *= layer_pdf
            if layer_pdf==0.: break

            
        # end iterating over layers
        #print ('{0:.5}'.format(path_throughput))
        amplitudes.append(path_throughput)
        pdfs.append(path_pdf)
        if abs(path_throughput) > EPSILON:
            amplitude += path_throughput
            non_zero_paths_ndx.append(path_ndx)

            if __DEBUG__:
                print (path_ndx, path)
                print (path_throughput, path_pdf)
                print ()

    # end iterating over paths
    print ()
    
    return amplitude, amplitudes, pdfs, non_zero_paths_ndx


    